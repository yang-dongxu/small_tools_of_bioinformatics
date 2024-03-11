split_txt() {
    local file_path="$1"
    local parts="$2"
    local output_dir="$3"
    local base_name=$(basename "$file_path")
    local stem_name="${base_name%.*}"

    local total_lines=$(wc -l < "$file_path")
    local header_lines=1
    local data_lines=$((total_lines - header_lines))
    local lines_per_part=$((data_lines / parts))
    local start_line=2
    local suf="${file_path##*.}"
    
    for i in $(seq 1 $parts); do
        if [ $i -eq $parts ]; then
            # For the last file, just copy to the end to avoid missing any lines due to rounding
            sed -n "1p;${start_line},\$p" "$file_path" > "${output_dir}/${stem_name}.part_${i}.${suf}"
        else
            local end_line=$((start_line + lines_per_part - 1))
            sed -n "1p;${start_line},${end_line}p" "$file_path" > "${output_dir}/${stem_name}.part_${i}.${suf}"
            start_line=$((end_line + 1))
        fi
    done
}

usage="$0 <input_txt> <parts> [odir]"
if [ $# -lt 2 ] || [ $# -gt 3 ] ;  then
    echo "Not enough or too much parameters input! "
    echo "usage:$usage"
    exit 1
fi
input_csv=$1
parts=$2

if [ $# -eq 3 ]; then
    odir=$3
else
    odir="./"
    echo "No odir specified! use `pwd` as default! "
fi

split_txt $input_csv $parts $odir
