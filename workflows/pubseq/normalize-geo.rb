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

# ---- For each metadata JSON file
Dir.new(GLOBAL.path).entries.select {|s| s =~/json$/}.each do |fn|
  next if fn == "state.json"
  jsonfn = GLOBAL.path+"/"+fn
  json = JSON.parse(File.read(jsonfn))
  # p json
  meta = OpenStruct.new(json)
  # ---- Step 1: find and normalize country
  collection_location = meta.sample['collection_location']
  # if text.match('countryname',collection_location)
  #   p
  Zlib::GzipReader.open('../../data/wikidata/countries.tsv.gz').each_line {|line|
    p line
  }
end
