#!/bin/bash
#
# This script converts .org files to .html so these generated
# files do not have to live in the git repo.

echo "Convert $1 from .org to .html"

guix environment --ad-hoc emacs-minimal emacs-org emacs-htmlize -- emacs -batch -visit $1 -eval "(progn (require 'org) (let ((org-export-htmlize-output-type 'css)) (org-html-export-to-html nil nil nil t nil)))"
