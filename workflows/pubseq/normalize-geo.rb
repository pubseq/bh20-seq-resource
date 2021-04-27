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

# =======================================================================
# ---- Actual processing starts here
# ---- For each metadata JSON file

# ---- Fetch state.json file
state = JSON.parse(File.read(GLOBAL.path+"/state.json"))

# Dir.new(GLOBAL.path).entries.select {|s| s =~/json$/}.each do |fn|
# next if fn == "state.json"

state.keys.each do |id|
  fn = id+".json"
  jsonfn = GLOBAL.path+"/"+fn
  json = JSON.parse(File.read(jsonfn))
  meta = OpenStruct.new(json)

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

  if not match_place.call(places_uris, lambda { |n| location =~ /:.*#{n}/ })
    match_place.call(places_uris, lambda { |n| address =~ /#{n}/ })
  end

  meta.warnings.push "Failed to normalize location" if not meta.sample['collection_location']

  # ---- Write new meta file
  json = JSON::pretty_generate(meta.to_h)
  print json,"\n" if GLOBAL.verbose
  File.write(GLOBAL.out+"/"+fn,json)
end
