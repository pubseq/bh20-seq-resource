A workflow to generate a phylogeny that can be visualized using [auspice](https://github.com/urbanslug/auspice).  
Expects a multi-fasta file path at [pggb_fasta][1] and generates a tree in `json` format.

#### Dependencies

Depends on:
 - [pggb](https://github.com/pangenome/pggb/blob/master/pggb)
   * [wfmash](https://github.com/ekg/wfmash)
   * [seqwish](https://github.com/ekg/seqwish)
   * [smoothxg](https://github.com/pangenome/smoothxg)
   * [odgi](https://github.com/vgteam/odgi)

 - [taxophages](https://github.com/urbanslug/taxophages/)
   * Clone and run with `python main.py ...`

 - [augur](https://github.com/nextstrain/augur)


#### Running

Expects that taxophages is cloned in a previous dir but you can update the path [main_py_script][2] to wherever it is.

Run the phylogeny workflow with the bleow after specifying your path to [pggb_fasta][1].

```bash
R_PACKAGES="${HOME}/RLibraries" \     # a directory holding R packages. Needed if R packages installed using install.packages on server e.g https://github.com/urbanslug/taxophages/blob/master/scripts/deps.R
TAXOPHAGES_ENV=server \               # helps taxophages figure out where it is being ran
AUGUR_RECURSION_LIMIT=30000 \         # augur isn't used to working with so many nested values
cwltool --preserve-entire-environment --no-container phylogeny.cwl clado-job.yml
```

Alternatively run any workflow with
```
cwltool --no-container <workflow>.cwl clado-job.yml
```

[1]: clado-job.yml#L8
[2]: clado-job.yml#L28
