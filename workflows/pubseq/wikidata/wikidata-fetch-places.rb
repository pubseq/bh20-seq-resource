#!/usr/bin/env ruby
#
# Use a SPARQL query to fetch Wikidata places
#
# Note: requires Ruby 3.x. Older Ruby gives a syntax error
#
# You may need to set
#
#   export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
#
# See also

raise "Currently not used"

require 'net/http'
require 'json'
require 'ostruct'
require 'erb'
require 'pp'

MAX=10

USER_AGENT =  {'User-Agent': 'genenetworkCrawler/1.0 (covid-19.genenetwork.org; pjotr.public821@thebird.nl) genenetworkCrawler/1.0', "Accept": "text/csv"}

SPARQL_HEADER="
prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
prefix dc: <http://purl.org/dc/terms/>
prefix schema: <https://schema.org/>
"

# Build a SPARQL query, submit and return results. Apply transform
# lambda when passed in
def sparql query, transform = nil

  api_url = 'https://query.wikidata.org/sparql'
  response = Net::HTTP.get_response(URI.parse(api_url),USER_AGENT)
  data = JSON.parse(response.body,symbolize_names: true)
  data => { head: { vars: }, results: { bindings: results} } # Ruby3 destructuring
  vars = vars.map { |v| v.to_sym }
  results.map { |rec|
    # return results after transforming to a Hash and applying the
    # optional transform lambda. Note the transform can not only
    # reduce results, or create an array, but also may transform into
    # an OpenStruct.
    res = {}
    vars.each { |name| res[name] = rec[name][:value] }
    if transform
      transform.call(res)
    else
      res
    end
  }
end

start = 0
num = MAX
begin
  query = "
SELECT DISTINCT ?place ?placename ?country ?coor ?population where {
  ?place wdt:P17 ?country ;
         wdt:P625 ?coor ;
         wdt:P1082 ?population .
  FILTER (?population > 9999)
  ?place rdfs:label ?placename .
  FILTER (lang(?placename)='en')
} LIMIT #{num} OFFSET #{start}
"
  list = sparql(query) # , lambda { |rec| rec[:id] })
  list.each do | l |
    print(l,"\n")
  end
  $stderr.print("#{start}-#{start+list.size}:#{list.first}\n") # show progress
  start += num
  exit 1
end while list.size == MAX
