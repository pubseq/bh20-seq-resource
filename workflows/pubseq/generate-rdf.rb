#! /usr/bin/env ruby
#
# -*- coding: UTF-8 -*-
#
# This script transforms pass2 JSON to JSON-LD (ready for RDF)
# See also https://github.com/pubseq/bh20-seq-resource/doc/blog/covid19-pubseq-update-rdf.org
#
# Author:: Pjotr Prins
# License:: MIT
#
# Copyright (C) 2021 Pjotr Prins <pjotr.prins@thebird.nl>
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

  o.on("--progress", "Show progress") do |p|
    options[:progress] = true
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
