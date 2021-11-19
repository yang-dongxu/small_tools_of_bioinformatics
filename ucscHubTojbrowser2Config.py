
## import librarise
import re
import os
import sys
import click
import requests
import json
import pathlib
import logging

# define Class Track
# track has format like below:
# track H3K9me3_M2
# type bigWig
# bigDataUrl datasets/public/H3K9me3_mouse_earlyembryo_2018_NCB/H3K9me3_Oocyte_M2.fc.bw
# shortLabel H3K9me3 ChIP-Seq in Oocyte
# longLabel H3K9me3 ChIP-Seq in Oocyte
# viewLimits 1:30
# color  152,251,152
# visibility full
class Track:
    def __init__(self,text):
        self.text = text
        self.track = ""
        self.type = ""
        self.bigDataUrl = ""
        self.shortLabel = ""
        self.longLabel = ""
        self.viewLimits = ""
        self.color = "0,0,255" ## defalut color is blue in rgb
        self.visibility = "full" ## defalut visibility is full
        self.parse()
    def parse(self):
        # each attribute is store in line start with corresponding key, and value is all the text after key
        # get attributes
        for line in self.text.split("\n"):
            if line.strip().startswith("track"):
                self.track = " ".join(line.split()[1:]).strip()
            elif line.strip().startswith("type"): 
                self.type = " ".join(line.split()[1:]).strip()
            elif line.strip().startswith("bigDataUrl"):
                self.bigDataUrl = " ".join(line.split()[1:]).strip()
            elif line.strip().startswith("shortLabel"):
                self.shortLabel = " ".join(line.split()[1:]).strip()
            elif line.strip().startswith("longLabel"):
                self.longLabel = " ".join(line.split()[1:]).strip()
            elif line.strip().startswith("viewLimits"):
                self.viewLimits = " ".join(line.split()[1:]).strip()
            elif line.strip().startswith("color"):
                self.color = " ".join(line.split()[1:]).strip()
            elif line.strip().startswith("visibility"):
                self.visibility = " ".join(line.split()[1:]).strip()
        return self

    # process to_dict for type bigWig
    def __to_dict_bigWig(self, remote = False, remote_url = ""):
        if remote:
            uri = "/".join(remote_url.split("/")[:-1] + [self.bigDataUrl])
        else:
            uri = self.bigDataUrl
        attrs = {}
        attrs["type"] = "QuantitativeTrack"
        attrs["trackId"] = self.track
        attrs["name"] = self.shortLabel
        attrs["adapter"] = {
            "type": "BigWigAdapter",
            "bigWigLocation":{
                "uri": uri,
                "locationType": "UriLocation"
            }
            
        }
        attrs["displays"] = [{}]
        dtype = "LinearWiggleDisplay"
        attrs["displays"][0]["type"] = dtype
        attrs["displays"][0]["displayId"] = f"{self.track}-{dtype}"
        if self.viewLimits:
            attrs["displays"][0]["minScore"] = float(self.viewLimits.split(":")[0])
            attrs["displays"][0]["maxScore"] = float(self.viewLimits.split(":")[1])
        if self.color:
            attrs["displays"][0]["renderers"] = {}
            render = attrs["displays"][0]["renderes"]
            render["XYPlotRenderer"] = {
                "type": "XYPlotRenderer",
                "posColor" : f"rgb({self.color})",
                "negColor" : f"rgb({self.color})"
            }
        ## if display is blank, drop it
        if not attrs["displays"][0]:
            del attrs["displays"]
        return attrs

    # process to_dict for type bigBed
    def __to_dict_bigBed(self, remote = False, remote_url = ""):
        if remote:
            uri = "/".join(remote_url.split("/")[:-1] + [self.bigDataUrl])
        else:
            uri = self.bigDataUrl
        attrs = {}
        attrs["type"] = "FeatureTrack"
        attrs["trackId"] = self.track
        attrs["name"] = self.shortLabel
        attrs["adapter"] = {
            "type": "BigBedAdapter",
            "bigBedLocation":{
                "uri": uri,
                "locationType": "UriLocation"
            }
        }
        attrs["displays"] = [{}]
        dtype = "LinearBasicDisplay"
        attrs["displays"][0]["type"] = dtype
        attrs["displays"][0]["displayId"] = f"{self.track}-{dtype}"
        if self.viewLimits:
            attrs["displays"][0]["minScore"] = float(self.viewLimits.split(":")[0])
            attrs["displays"][0]["maxScore"] = float(self.viewLimits.split(":")[1])
        if self.color:
            attrs["displays"][0]["color"] = f"rgb({self.color})"
        ## if display is blank, drop it
        if not attrs["displays"][0]:
            del attrs["displays"]
        return attrs

    # fit attributes into a dict
    def to_dict(self, remote = False, remote_url = ""):
        if self.type == "bigWig":
            attrs = self.__to_dict_bigWig(remote, remote_url)
        elif self.type == "bigBed":
            attrs = self.__to_dict_bigBed(remote, remote_url)
        else:
            raise Exception(f"type {self.type} is not supported")
        
        return attrs



# defind a class contain ucscHub information
class UcscHub:
    def __init__(self, url) -> None:
        self.url = url
        if self.url.startswith('http://') or self.url.startswith('https://'):
            self.type = "remote"
        elif self.url.startswith('file://'):
            self.type = "local"
        else:
            raise ValueError("url should start with http:// or https:// or file://")

        self.hub_text = ""
        self.hub_name = ""
        self.short_label = ""
        self.long_label = ""

        self.genomeFile = ""
        self.genome_text = ""
        self.genome = ""

        self.trackDb = ""
        self.tracks_text = ""
        self.tracks = []

        self.get_hub()
        self.get_genome(self.get_hub())
        self.get_tracks()
        self.parse_tracks()
        return 

    ## get hub content, and parse it, and return genomeFile name
    def get_hub(self):
        # if remote, download it
        if self.type == "remote":
            r = requests.get(self.url)
            self.hub = r.text
        # if local, read it
        elif self.type == "local":
            with open(self.url[7:], 'r') as f:
                self.hub = f.read()

        ## hubname: first line
        self.hub_name = self.hub.splitlines()[0].split()[1]
        ## shortlabel: line start with shortLabel
        self.short_label = re.findall(r'\bshortLabel\s+(.*)\b', self.hub)[0]
        ## longlabel: line start with longLabel
        self.long_label = re.findall(r'\blongLabel\s+(.*)\b', self.hub)[0]
        ## genome: line start with genomesFile
        genomeFile = re.findall(r'\bgenomesFile\s+(.*)\b', self.hub)[0]
        return genomeFile

    ## get genomeFile content, and parse it, and trackDB name
    def get_genome(self, genomeFile):
        ## replace url last term by genomeFile name
        self.genomeFile = "/".join(self.url.split('/')[:-1] + [genomeFile])
        ## if remote, download it
        if self.type == "remote":
            r = requests.get(self.genomeFile)
            self.genome_text = r.text
        # if local, read it
        elif self.type == "local":
            with open(self.genomeFile, 'r') as f:
                self.genome_text = f.read()
        # genome: first line
        self.genome = self.genome_text.splitlines()[0].split()[1]
        # tracks: line start with trackDb
        trackDb = self.genome_text.splitlines()[1].split()[1]
        # replace url last term by trackDb name
        self.trackDb = "/".join(self.url.split('/')[:-1] + [trackDb])

        return None

    ## get trackDb content, and parse it, and return tracks
    def get_tracks(self):
        # if remote, download it
        if self.type == "remote":
            r = requests.get(self.trackDb)
            self.tracks_text = r.text
        # if local, read it
        elif self.type == "local":
            with open(self.trackDb, 'r') as f:
                self.tracks_text = f.read()
        return self

    ## parse tracks_text, and return tracks
    def parse_tracks(self):
        # each track is a part of text, the first line start with "track", and the last line end is followed by a blank line
        # split trakcs_text by re
        tracks_text = self.tracks_text.split("track")
        # remove the first element, which is empty
        tracks_text.pop(0)
        # for each element, add start with "track" to the front, and remove blanks at the end
        tracks_text = ["track" + track for track in tracks_text]
        tracks_text = [track.rstrip() for track in tracks_text]
        # for each element, create a Track object
        self.tracks = [Track(track) for track in tracks_text]
        return self.tracks

    ## return a list, contain all tracks dict, with modify genome name
    def to_dict(self, remote = False):
        tracks = []
        for track in self.tracks:
            new_track = track.to_dict(remote, self.url)
            new_track["assemblyNames"] = [self.genome]
            new_track["category"] = [self.hub_name, self.long_label]
            # del new_track["displays"]
            tracks.append(new_track)
        return tracks




# define a function to accept params:
# -a --ucsc_hub_url, can be a list, or multi terms separated by comma; if it's a remote url, start with http://; or local file, start with file:// \n
# -b --base_jbrowse_config_url; if it's a remote url, start with http://; or local file, start with file:// \n
# -r --remote, whether save output file in remote url format or simply save it as described in ucscHub file 
# -o --output, output file name, default is stdout
@click.option("-a","--ucsc_hub_url",type=str, multiple=True,required=True, help="ucsc hub url, tsv format \n")
@click.option("-b","--base_jbrowse_config_url",type=str,required=True, help="base jbrowse config url, json format \n")
@click.option("-r","--remote", is_flag=True, help="remote url, json format \n")
@click.option("-o","--output",type=click.File("w"),required=False, default=sys.stdout, help="output file name \n")
@click.command()
def process(ucsc_hub_url,base_jbrowse_config_url,remote,output):
    ## tidy ucsc_hub_url
    ucsc_hub_url = list(set([j for i in ucsc_hub_url for j in i.split(",")]))

    ## create a list of UcscHub objects
    ucscHubs = [UcscHub(url) for url in ucsc_hub_url]


    ## if base_jbrowse_config_url is remote, download it
    if base_jbrowse_config_url.startswith('http://') or base_jbrowse_config_url.startswith('https://'):
        r = requests.get(base_jbrowse_config_url)
        baseJbrowseConfig = r.text
    # if local, read it
    elif base_jbrowse_config_url.startswith('file://'):
        with open(base_jbrowse_config_url[7:], 'r') as f:
            baseJbrowseConfig = f.read()
    baseJbrowseConfig = json.loads(baseJbrowseConfig)

    ## append ucscHubs to baseJbrowseConfig["tracks"]
    for ucscHub in ucscHubs:
        baseJbrowseConfig["tracks"].extend(ucscHub.to_dict(remote))
    
    ## output baseJbrowseConfig
    json.dump(baseJbrowseConfig, output, indent=4)

    return ucscHubs, baseJbrowseConfig, remote, output



if __name__ == "__main__":
    process()