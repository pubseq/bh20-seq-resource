import sys

reads = $(inputs.readsFA)
subjects = $(inputs.subjects)

for i, r in enumerate(reads):
    with open(r["path"], "rt") as fa:
        fa.readline()
        print(">"+subjects[i])
        data = fa.read(8096)
        while data:
            sys.stdout.write(data)
            data = fa.read(8096)
