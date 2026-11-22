# Animate the dynamics of a food web

Writes a GIF in which the trajectories from
[`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md)
unfold over time, one frame per time step shown, with the panels and
axes of
[`plot.food_web()`](https://provetelab.org/comdynamite/reference/plot.food_web.md).
It needs the package gifski.

## Usage

``` r
animate_food_web(
  x,
  file = tempfile(fileext = ".gif"),
  frames = 100,
  delay = 0.08,
  width = 640,
  height = 640,
  ...
)
```

## Arguments

- x:

  Result of
  [`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md).

- file:

  Path of the GIF to write, one character string; by default a temporary
  file, whose path is returned.

- frames:

  Number of frames, taken at evenly spaced rows of `x`: one whole
  number, 2 or more, that cannot exceed the number of rows.

- delay:

  Seconds between frames, one number greater than 0.

- width, height:

  Size of the GIF in pixels, each one whole number, 1 or more.

- ...:

  Other arguments of
  [`plot.food_web()`](https://provetelab.org/comdynamite/reference/plot.food_web.md),
  such as `panels` or `colours`.

## Value

The path of the GIF, invisibly.

## See also

[`plot.food_web()`](https://provetelab.org/comdynamite/reference/plot.food_web.md)

Other coevolution in food webs:
[`food_web_dynamics()`](https://provetelab.org/comdynamite/reference/food_web_dynamics.md),
[`food_web_params()`](https://provetelab.org/comdynamite/reference/food_web_params.md),
[`plot.food_web()`](https://provetelab.org/comdynamite/reference/plot.food_web.md)

## Examples

``` r
# \donttest{
if (requireNamespace("gifski", quietly = TRUE)) {
  x <- food_web_dynamics(seq(0, 3000, by = 30), dt = 0.5)
  gif <- animate_food_web(x, tempfile(fileext = ".gif"), frames = 30)
  file.exists(gif)
}
#> food_web_dynamics() implements the method of:
#>   McPeek, M. A. (2017). Evolutionary Community Ecology. Princeton University Press. https://doi.org/10.23943/princeton/9780691088778.001.0001
#>   McPeek, M. A. (2017). The ecological dynamics of natural selection: traits and the coevolution of community structure. The American Naturalist, 189, E91-E117. https://doi.org/10.1086/691101
#> Please cite these papers together with comdynamite (citation("comdynamite")).
#> This message is shown once per session; options(comdynamite.cite = FALSE) turns it off.
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> `geom_line()`: Each group consists of only one observation.
#> ℹ Do you need to adjust the group aesthetic?
#> [1] TRUE
# }
```
