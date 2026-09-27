# Synthetic reverse-primer regression fixture

These sequences test software behavior and do not describe a biological assay.
The primer is supplied as the 40-base oligo in 5′ → 3′ direction. Its literal
reverse complement is `TAGGTCGAACTGCTAGCATCGATCGCTAAGCTTGCAACGT`.

Each binding region has 40 C bases before it and 40 G bases after it. Every
variant is included in forward and reverse-complemented orientation. The
primer role remains reverse in both cases, so tests also check that role
metadata does not restrict which FASTA strand is searched.

| Variant | Expected differences in primer orientation |
| --- | --- |
| perfect | None |
| first | 1:A>T (5′ tip) |
| last | 40:A>C (3′ tip) |
| both | 1:A>T,40:A>C; BLAST clips both ends |
| internal | 16:C>A |
| insertion | One inserted G, anchored after position 18 |
| deletion | 19:T>- |

The inserted G differs from its adjacent bases A and T to avoid equivalent
gap placements within a repeated-base run.
