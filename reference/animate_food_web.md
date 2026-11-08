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
  file = "food_web.gif",
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

  Path of the GIF to write.

- frames:

  Number of frames, taken at evenly spaced rows of `x`; it cannot exceed
  the number of rows.

- delay:

  Seconds between frames.

- width, height:

  Size of the GIF in pixels.

- ...:

  Other arguments of
  [`plot.food_web()`](https://provetelab.org/comdynamite/reference/plot.food_web.md),
  such as `panels` or `col`.

## Value

The path of the GIF, invisibly.

## See also

[`plot.food_web()`](https://provetelab.org/comdynamite/reference/plot.food_web.md)

## Examples

``` r
# \donttest{
if (requireNamespace("gifski", quietly = TRUE)) {
  x <- food_web_dynamics(seq(0, 3000, by = 30), dt = 0.5)
  gif <- animate_food_web(x, tempfile(fileext = ".gif"), frames = 30)
  file.exists(gif)
}
#> [1] TRUE
# }
```
