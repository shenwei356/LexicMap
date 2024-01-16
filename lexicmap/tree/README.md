## tree

This package is modified from https://github.com/armon/go-radix, for support the bit-packed k-mers.
Some new methods are added to support querying keys shared with a prefix, not the original `LongestPrefix`.

This package requires indexed and query k-mers having the same K value.
