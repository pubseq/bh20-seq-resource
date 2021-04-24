#! /usr/bin/env ruby
#
# -*- coding: UTF-8 -*-
#
# Metadata validation routine - does one JSON file. At this stage this
# is mostly a debugging tool
#
# Author:: Pjotr Prins
# License:: MIT
#
# Copyright (C) 2021 Pjotr Prins <pjotr.prins@thebird.nl>
#
#

TOOL=File.basename($0)

GEMPATH = File.dirname(__FILE__) + '/../../lib/ruby'
$: << File.join(GEMPATH,'lib/ruby/pubseq')

VERSION_FILENAME=File.join(GEMPATH,'VERSION')
VERSION = File.new(VERSION_FILENAME).read.chomp

require 'colorize'
require 'optparse'
require 'ostruct'
require 'fileutils'
require 'json'
require 'zlib'

options = { show_help: false, source: 'https://github.com/pubseq', version: VERSION+' (Pjotr Prins)', date: Time.now.to_s }

opts = OptionParser.new do |o|
  o.banner = "Usage: #{TOOL} [options] path"

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
  print "\nExample: ruby validate.rb MT810507.json -q|jq\n"
  exit 1
end

if RUBY_VERSION =~ /^[12]/
  $stderr.print "WARNING: #{TOOL} may not run properly on Ruby <3.x\n"
end

$stderr.print "Options: ",options,"\n" if !options[:quiet]

GLOBAL = OpenStruct.new(options)
$has_error = false


for fn in ARGV
  next if fn == "state.json"
  json = JSON.parse(File.read(fn))
  meta = OpenStruct.new(json)
  sample = OpenStruct.new(meta.sample)

  error = lambda { |msg|
    print(json.to_json,"\n")
    $stderr.print "ERROR: ".red,msg.red,"\n"
    $has_error = true
  }

  # ---- Check for location
  location = meta.sample['collection_location']
  error.call "Missing collection_location" if not location
  error.call "Collection_location <#{location}> not normalized" if location !~ /^http:\/\/www.wikidata.org\/entity\/Q/

  # ---- Dates
  error.call "Sample collection_date <#{sample.collection_date}> malformed" if sample.collection_date !~ /\d\d\d\d-\d\d-\d\d/

end

exit 1 if $has_error
