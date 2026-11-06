# Daphniids of the northeastern United States

Incidence of 13 daphniids in lakes of the EMAP Northeastern Lakes Survey
of 1991 to 1994 (first visit to each lake), rebuilt to follow Leibold et
al. (2010). The tables are public data of the US Environmental
Protection Agency, and the script is `data-raw/daphniids.R`.

## Usage

``` r
daphniids
```

## Format

List with `comm`, a matrix of lakes by species; `env`, a matrix with
colour, conductivity, dissolved organic carbon, total nitrogen, pH,
total phosphorus, mean depth, elevation and the ratio of nitrogen to
phosphorus; `biogeography`, a matrix of contrasts with values of 1, -1
and 0; and `tree`, a `phylo` object.

## Source

<https://archive.epa.gov/emap/archive-emap/web/html/nelakes.html>

## Details

The article gives 327 lakes and this rebuilt set has 314, of which 239
have some species. The tree is the topology of Figure 3a of the article,
without branch lengths. The survey lists *D. longiremis* and not the *D.
longispina* of the figure, which we took to be the same taxon. The
article does not give the area of each lake, so rivers in `biogeography`
are approximated by longitude and these contrasts differ from the
original ones: S1, east and west of the Hudson; S2, Adirondacks and
other lakes to the west of the Hudson; S3, east and west of the
Connecticut River, east of the Hudson; S4, Maine and other lakes east of
the Hudson; S5, New Hampshire and lakes to its west but east of the
Hudson; S6, New Hampshire and lakes in Rhode Island and Massachusetts.

## References

Leibold, M. A., Economo, E. P., & Peres-Neto, P. (2010). Metacommunity
phylogenetics: separating the roles of environmental filters and
historical biogeography. *Ecology Letters*, 13, 1290–1299.
[doi:10.1111/j.1461-0248.2010.01523.x](https://doi.org/10.1111/j.1461-0248.2010.01523.x)
