#! /usr/bin/env ruby
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
  o.on('--dir path',String, 'Path to JSON files') do |path|
    options[:path] = path
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

NORMALIZE_GEO_BANNER = "#{TOOL} #{VERSION} (Ruby #{RUBY_VERSION}) by Pjotr Prins 2021\n"
$stderr.print NORMALIZE_GEO_BANNER if !options[:quiet]

if options[:show_help]
  print opts
  exit 1
end

if RUBY_VERSION =~ /^[12]/
  $stderr.print "WARNING: #{TOOL} may not run properly on Ruby <3.x\n"
end

$stderr.print "Options: ",options,"\n" if !options[:quiet]

GLOBAL = OpenStruct.new(options)

# ---- So far it is boiler plate to set the environment
country_uri = {}
Zlib::GzipReader.open('../../data/wikidata/countries.tsv.gz',:encoding => 'UTF-8').each_line {|line|
  place,country,countryname,continent = line.split(/\t/)
  country_uri[countryname.strip] = place
}
country_alias_uri = {}
Zlib::GzipReader.open('../../data/wikidata/country_aliases.tsv.gz',:encoding => 'UTF-8').each_line {|line|
  place,country,countryname,aliases = line.split(/\t/)
  country_alias_uri[aliases.upcase.strip] = place
}

# ---- For each metadata JSON file
Dir.new(GLOBAL.path).entries.select {|s| s =~/json$/}.each do |fn|
  next if fn == "state.json"
  jsonfn = GLOBAL.path+"/"+fn
  json = JSON.parse(File.read(jsonfn))
  # p json
  meta = OpenStruct.new(json)
  match = lambda { |hash,location|
    loc = location.upcase.strip
    hash.each { |c,uri|
      if loc =~ /\W#{c}\W/
        meta.sample['collection_location'] = uri
        meta.sample['original_collection_location'] = location
        meta.sample['country'] = c
        meta.sample['wd:country'] = uri
        return true
      end
    }
    false
  }
  # ---- Step 1: find and normalize by country
  collection_location = meta.sample['collection_location']
  if collection_location
    loc = collection_location.strip
    if not match.call(country_uri,loc)
      # ---- Step 2: find and normalize by alias
      if not match.call(country_alias_uri,loc)
        # ---- Step 3: find and normalize using submitter address
        submitter_address = meta.submitter['submitter_address']
        if submitter_address
          match.call(country_uri,submitter_address)
        end
      end
    end
  end
  p meta
end
