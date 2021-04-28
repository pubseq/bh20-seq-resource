#! /usr/bin/env ruby
#
# -*- coding: UTF-8 -*-
#
# PubSeq normalize geo info using wikidata URIs. Based on earlier
# PubSeq scripts (written in Python)
#
# Author:: Pjotr Prins
# License:: MIT
#
# Copyright (C) 2021 Pjotr Prins <pjotr.prins@thebird.nl>
#
# This script takes a metadata file and transforms the geo/location
# fields into a normalized wikidata URI. This is described in
# https://github.com/pubseq/bh20-seq-resource/blob/master/doc/blog/covid19-pubseq-location-data.org
#

TOOL=File.basename($0)

GEMPATH = File.dirname(__FILE__) + '/../../lib/ruby'
$: << File.join(GEMPATH,'lib/ruby/pubseq')

VERSION_FILENAME=File.join(GEMPATH,'VERSION')
VERSION = File.new(VERSION_FILENAME).read.chomp

require 'optparse'
require 'ostruct'
require 'fileutils'
require 'json'
require 'zlib'

options = { show_help: false, source: 'https://github.com/pubseq', version: VERSION+' (Pjotr Prins)', date: Time.now.to_s }

opts = OptionParser.new do |o|
  o.banner = "Usage: #{TOOL} [options] path"
  o.on('--dir path',String, 'Path to JSON files [REQUIRED]') do |path|
    options[:path] = path
  end
  o.on('--out path',String, 'Dir to write to [REQUIRED]') do |path|
    options[:out] = path
  end

  o.separator ""

  o.on("-q", "--quiet", "Run quietly") do |q|
    # Bio::Log::CLI.trace('error')
    options[:quiet] = true
  end

  o.on("-v", "--verbose", "Run verbosely") do |v|
    options[:verbose] = true
  end

  o.on("--debug", "Show debug messages and keep intermediate output") do |v|
    # Bio::Log::CLI.trace('debug')
    options[:debug] = true
  end

  o.separator ""
  o.on_tail('-h', '--help', 'display this help and exit') do
    options[:show_help] = true
  end
end

opts.parse!(ARGV)

BANNER = "#{TOOL} #{VERSION} (Ruby #{RUBY_VERSION}) by Pjotr Prins 2021\n"
$stderr.print BANNER if !options[:quiet]

if options[:show_help]
  print opts
  exit 1
end

if RUBY_VERSION =~ /^[12]/
  $stderr.print "WARNING: #{TOOL} may not run properly on Ruby <3.x\n"
end

$stderr.print "Options: ",options,"\n" if !options[:quiet]

GLOBAL = OpenStruct.new(options)

raise "--dir directory is required" if not GLOBAL.path
raise "--out directory is required" if not GLOBAL.out


# =======================================================================
# ---- So far it is boiler plate to set the environment
# ---- continue reading the tables from Wikidata

country_uri = [] # tuples
Zlib::GzipReader.open('../../data/wikidata/countries.tsv.gz',:encoding => 'UTF-8').each_line {|line|
  place,country,countryname,continent = line.split(/\t/)
  country_uri.push [countryname.upcase.strip,place]
}

explicit = [
  ["USA","http://www.wikidata.org/entity/Q30"],
]

country_uris = country_uri.sort_by { |t| t[0].length }.reverse + explicit

alias_uri = []
Zlib::GzipReader.open('../../data/wikidata/country_aliases.tsv.gz',:encoding => 'UTF-8').each_line {|line|
  place,country,countryname,aliases = line.split(/\t/)
  alias_uri.push [aliases.upcase.strip,place]
}
alias_uris = alias_uri.sort_by { |t| t[0].length }.reverse

places_uri = []
Zlib::GzipReader.open('../../data/wikidata/places.csv.gz',:encoding => 'UTF-8').each_line {|line|
  placename,place,country,coor,population = line.split(/,/)
  places_uri.push [placename.strip,place,country,coor,population]
}
Zlib::GzipReader.open('../../data/wikidata/regions.csv.gz',:encoding => 'UTF-8').each_line {|line|
  placename,place,country,coor,population = line.split(/,/)
  places_uri.push [placename.strip,place,country,coor,population]
}
places_uris = places_uri.sort_by { |t| t[0].length }.reverse

# ========================================================================
# ---- Actual processing starts here
# ---- For each metadata JSON file

# ---- Fetch state.json file
state = JSON.parse(File.read(GLOBAL.path+"/state.json"))

# Dir.new(GLOBAL.path).entries.select {|s| s =~/json$/}.each do |fn|
# next if fn == "state.json"

state.keys.each do |id|
  next if not state[id]['valid']
  fn = id+".json"
  jsonfn = GLOBAL.path+"/"+fn
  json = JSON.parse(File.read(jsonfn))
  meta = OpenStruct.new(json)

  warn = lambda { |msg|
    msg += " for #{id}"
    meta.warnings.push msg
    $stderr.print "WARNING: #{msg}\n" if GLOBAL.verbose or GLOBAL.debug
  }

  # ---- GEO Normalisation -----------------------------------------------
  # ==== Normalize by country using uploader location fields for
  #      collection_location and submitter_address

  # ---- Check for location
  location = meta.sample['collection_location']
  # meta.warnings.push "Missing collection_location" if not location - caught in yamlfa

  # ---- Update orginal_collection_location
  if not meta.sample['original_collection_location']
    meta.sample['original_collection_location'] = location if location
    meta.sample.delete('collection_location')
  end

  address = meta.submitter['submitter_address']

  match = lambda { |wds,find|
    wds.each { |name,uri|
      if find.call(name)
        meta.sample['collection_location'] = uri
        meta.sample['country'] = name
        meta.sample['wd:country'] = uri
        return true
      end
    }
    false
  }

  location = ( location ? location.upcase.strip : "UNDEF")
  address = ( address ? address.upcase.strip : "UNDEF")
  # GenBank countries are formed as "USA: MD"
  if not match.call(country_uris, lambda { |n| location =~ /^#{n}:/ui or location.end_with?(n) })
    # Try aliases
    if not match.call(alias_uris, lambda { |n| location.start_with?(n) })
      # Try submitter address
      if not match.call(country_uris, lambda { |n| address =~ /#{n}$/ui })
        # Try aliases on submitter address
        match.call(alias_uris, lambda { |n| address.end_with?(n) })
      end
    end
  end

  # ---- Refine location using the places/regions information
  wd_country = meta.sample['wd:country']
  countryname = meta.sample['country']

  match_place = lambda { |places,find|
    places.each { |name,uri,country_uri,coor,population|
      if wd_country == country_uri and find.call(name)
        meta.sample['collection_location'] = uri
        meta.sample['place'] = name
        return true
      end
    }
    false
  }

  location = meta.sample['original_collection_location']
  address = meta.submitter['submitter_address']
  strain = meta.virus['virus_strain']
  if not match_place.call(places_uris, lambda { |n| location =~ /:.*#{n}/ })
    if not match_place.call(places_uris, lambda { |n| address =~ /#{n}/ })
      match_place.call(places_uris, lambda { |n| strain =~ /#{n}/ })
    end
  end
  meta.sample.delete('wd:country') # implicit in Wikidata location

  warn.call "Failed to normalize location" if not meta.sample['collection_location']

  # ---- Normalise sequencing method -------------------------------------
  known_sequencers = [
    [ "HiSeq.?1000", "http://www.ebi.ac.uk/efo/EFO_0004204" ],
    [ "HiSeq.?2000", "http://www.ebi.ac.uk/efo/EFO_0004203" ],
    [ "HiSeq.?2500", "http://www.ebi.ac.uk/efo/EFO_0008565" ],
    [ "HiSeq.?3000", "http://www.ebi.ac.uk/efo/EFO_0008564" ],
    [ "HiSeq.?4000", "http://www.ebi.ac.uk/efo/EFO_0008563" ],
    [ "HiSeq.?100", "http://www.ebi.ac.uk/efo/EFO_0008635" ],
    [ "HiSeq.?X", "http://www.ebi.ac.uk/efo/EFO_0008567" ],
    [ "NextSeq.?500", "http://www.ebi.ac.uk/efo/EFO_0009173" ],
    [ "NextSeq", "http://www.ebi.ac.uk/efo/EFO_0008566" ],
    [ "MiSeq", "http://www.ebi.ac.uk/efo/EFO_0004205" ],
    [ "MiniSeq", "http://www.ebi.ac.uk/efo/EFO_0008636" ],
    [ "NovaSeq", "http://www.ebi.ac.uk/efo/EFO_0008637" ],
    [ "Genome.Analyzer.IIx", "http://www.ebi.ac.uk/efo/EFO_0004202" ],
    [ "Genome.Analyzer.II", "http://www.ebi.ac.uk/efo/EFO_0004201" ],
    [ "Genome.Analyzer", "http://www.ebi.ac.uk/efo/EFO_0004200" ],
    [ "Illumina", "http://purl.obolibrary.org/obo/OBI_0000759" ],
    [ "MinION", "http://www.ebi.ac.uk/efo/EFO_0008632" ],
    [ "Nanopore", "http://purl.obolibrary.org/obo/NCIT_C146818" ],
    [ "GridION", "http://www.ebi.ac.uk/efo/EFO_0008633" ],
    [ "PremethION", "http://www.ebi.ac.uk/efo/EFO_0008634" ],
    [ "Oxford.?Nanopore", "http://purl.obolibrary.org/obo/NCIT_C146818" ],
    [ "Ion.?Torrent", "http://purl.obolibrary.org/obo/NCIT_C125894" ],
    [ "Thermo.?Fisher", "http://purl.obolibrary.org/obo/NCIT_C125894" ],
    [ "Sanger", "http://purl.obolibrary.org/obo/NCIT_C19641" ],
    [ "MGISEQ.?2000", "http://virtual-bh/MGISEQ2000" ],
    [ "PacBio.?RS", "http://www.ebi.ac.uk/efo/EFO_0008631" ],
    [ "PacBio", "http://www.ebi.ac.uk/efo/EFO_0008630" ],
    [ "454", "http://www.ebi.ac.uk/efo/EFO_0004431" ],
    [ "SOLiD", "http://www.ebi.ac.uk/efo/EFO_0004435" ],
  ]

  if meta.technology['sample_sequencing_technology']
    meta.technology['sample_sequencing_technology'].map! { |seq|
      index = known_sequencers.index { |pair|
        name,uri = pair
        ( seq =~ /#{name}/i ? uri : false )
      }
      if index
        known_sequencers[index][1]
      else
        warn.call "Failed to normalize sequencer <#{seq}>"
      end
    }
  end

  # ---- Normalise assembly method -------------------------------------
  known_assembly = [
    [ "spades", "http://purl.obolibrary.org/obo/GENEPIO_0001628" ], # de novo
    [ "de novo", "http://purl.obolibrary.org/obo/GENEPIO_0001628" ], # de novo
    [ "assembly", "http://purl.obolibrary.org/obo/GENEPIO_0001628" ], # de novo
    [ "bwa", "http://purl.obolibrary.org/obo/GENEPIO_0002028" ], # default mapped
  ]

  ap = meta.technology['alignment_protocol']
  if ap
    known_assembly.each do |pair|
      name,uri = pair
      if ap =~ /#{name}/
        meta.technology['assembly_method'] = uri
        break
      end
    end
  else
    warn.call "Missing alignment protocol"
  end
  am = meta.technology['assembly_method']
  if not am or am !~ /^http/
    meta.technology['assembly_method'] = "http://purl.obolibrary.org/obo/GENEPIO_0002028"
  end


  # ---- Wrap up ---------------------------------------------------------
  # ---- Write new meta file
  json = JSON::pretty_generate(meta.to_h)
  print json,"\n" if GLOBAL.verbose
  File.write(GLOBAL.out+"/"+fn,json)

  # ---- Update state file
  state[id]["warnings"] = meta.warnings
end

# ---- Write the state file
state_json = JSON::pretty_generate(state)
File.write(GLOBAL.out+"/state.json",state_json)
