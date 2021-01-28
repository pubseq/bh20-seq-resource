# Normalization steps

This library contains generic logic to normalize (string) data and
transforms strings to URIs.  It should be applicable to data from
any source (GenBank, ENA etc).

Important: missing data should be missing or None! Do not fill
in data by 'guessing'.

When data is malformed a warning should be logged and added to the
warning list. Functions should be small enough to return only 1
warning!

Pjotr Prins (c) 2021
