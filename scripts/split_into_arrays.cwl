cwlVersion: v1.1
class: ExpressionTool
requirements:
  InlineJavascriptRequirement: {}
inputs:
  dir:
    type: Directory
    loadListing: shallow_listing
outputs:
  fasta: File[]
  metadata: File[]
expression: |
  ${
  var dir = inputs.dir;
  var fasta = [];
  var metadata = [];
  dir.listing.sort(function(a, b) { return a.basename < b.basename; });
  for (var i = 0; i < dir.listing.length; i++) {
    if (dir.listing[i].basename.substr(-6) == ".fasta") {
      fasta.push(dir.listing[i]);
    }
    if (dir.listing[i].basename.substr(-5) == ".yaml") {
      metadata.push(dir.listing[i]);
    }
  }
  if (fasta.length != metadata.length) {
    throw "They dont match";
  }
  return {"fasta": fasta, "metadata": metadata};
  }
