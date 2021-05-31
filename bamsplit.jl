using XAM;
using ArgParse;
using DataStructures;
using JSON;
using BioAlignments;
using BGZFStreams;

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--input", "-i"
            help = "input bam file"
            dest_name = "ibam"
            required = true 
        "--snp", "-s"
            help = "input snp file, in vcf format"
            dest_name = "isnp"
            required = true
        "--oname1", "-a"
            help = "output name1 "
            dest_name = "obam1"
            default = "strain1.bam"
        "--oname2", "-b"
            help = "output name2"
            dest_name = "obam2"
            default = "strain2.bam"
        "--mapq", "-m"
            help = "mapq threod "
            dest_name = "mapq"
            default = 30
            arg_type = Int
        "--quality","-q"
            help="quality of read base"
            dest_name="quality"
            default = 30
            arg_type = Int
    end
    return parse_args(s)
end

function get_snp(snp::String) 
    if ! isfile(snp)
        println("File not exist! $(snp)")
        exit()
    end
    p=DefaultDict{Int,Tuple{String,String} }(("*","*"))
    snp_stat=DefaultDict{String, DefaultDict{Int,Tuple{String,String} } }( p )
    f = open(snp)
    # line_number
    # read till end of file
    while ! eof(f)  
        line=readline(f)
        if startswith(line,"#")
            continue
        end
        # read a new / next line for every iteration           
        s = split(line)
        if startswith(s[1], "chr")
            chr=s[1]
        else
            chr="chr$(s[1])" 
        end   
        snp_stat[chr][parse(Int,s[2])]=(s[4],s[5])
    end
    close(f)
    return snp_stat
end

function quality_good(q,quality)
    p=Int(q) ## XAM 已经自动识别了phred33+， 不需要自己再手动变换
    if p>quality
        return true 
    else
        return false
    end
end

function process_record!(record::XAM.BAM.Record,stats::AbstractDict,snp_stat::AbstractDict,mapq=30,quality=30)
    mapq_i=BAM.mappingquality(record)
    if mapq_i < mapq
        return
    end
    align=BAM.alignment(record)
    seq=BAM.sequence(record)
    qs=BAM.quality(record)
    chr=BAM.refname(record)
    name=BAM.tempname(record)
    if ! haskey(stats,name)
        stats[name]=Dict{BAM.Record, Vector{Int}}(record=>[0,0,0])
    else
        stats[name][record]=[0,0,0]
    end
    #stats[name][record]=[0,0,0]
    for i in 1:length(seq)
        pos, cg= BioAlignments.seq2ref(align,i)
        base="$(seq[i])"
        q=qs[i]
        if ! quality_good(q,quality)
            continue
        end
        change=snp_stat[chr][pos]
        change[1]::String
        base::String
        if change[1] == base
            stats[name][record][1]+=1
        elseif change[2] == base
            stats[name][record][2]+=1
        else
            stats[name][record][3]+=1
        end
    end
end

function output(obam1,obam2,stats,header)
    fo1 = BAM.Writer(BGZFStream(open(obam1, "w"), "w"),header)
    fo2 = BAM.Writer(BGZFStream(open(obam2, "w"), "w"),header)

    function output_(f, records::Vector{XAM.BAM.Record})
        for record in records
            write(f,record)
        end
    end


    for (name, data) in stats
        vote=[0,0,0]::Vector{Int}
        records= Vector{BAM.Record}()
        for record in collect(keys(data))
            push!(records,record)
            vote += data[record]
        end
        all_vote=vote[1]+vote[2]
        if all_vote <2  ## snp 数量太少，可能是噪音，排除掉
            #println("No snp read: $(name)")
            continue
        elseif vote[1]/all_vote > 0.66 ##strain1
            output_(fo1,records)
        elseif vote[2]/all_vote > 0.66 ##strain2
            output_(fo2,records)
        else ## mixed
            println("mixed read: $(name)")
            println(vote)
            continue
        end
    end
    close(fo1)
    close(fo2)
end

function bam_split(ibam,snp,obam1,obam2,mapq=30,quality=30)
    snp_stat=get_snp(snp)
    ibam=open(BAM.Reader, ibam)
    header=XAM.BAM.header(ibam)

    #tmp=DefaultDict{XAM.BAM.Record,Vector{Int} }([0,0,0]) ## strain1, strain2, other
    #stats=DefaultDict{String,DefaultDict{XAM.BAM.Record,Vector{Int} } }( DefaultDict{XAM.BAM.Record,Vector{Int}}() )
    stats=Dict{String, Dict{BAM.Record, Vector{Int}} }()

    #record=XAM.BAM.Record()
    #while !eof(ibam)
    for record in ibam
        #empty!(record)
        #read!(ibam, record)
        process_record!(record,stats,snp_stat,mapq,quality)
    end
    output(obam1,obam2,stats,header)
end

function main()
    args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end
    bam_split(args["ibam"],args["isnp"],args["obam1"],args["obam2"],args["mapq"],args["quality"])
end

main()
