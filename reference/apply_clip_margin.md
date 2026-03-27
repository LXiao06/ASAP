# Apply Time Margin to Clip Rows and Validate Against WAV Duration

Expands `start_time` and `end_time` in a clip data frame by the
requested margin, then validates each clip against its source WAV file
duration (read via header only). Clips whose adjusted `start_time < 0`
or adjusted `end_time > wav_duration` are silently dropped and counted.

## Usage

``` r
apply_clip_margin(clips, wav_dir, margin)
```

## Arguments

- clips:

  Prepared clip data frame (must have `filename`, `start_time`,
  `end_time`, and optionally `day_post_hatch`).

- wav_dir:

  Root directory containing source WAV files.

- margin:

  Numeric vector of length 2: seconds to prepend before start and append
  after end. Both values must be non-negative.

## Value

A list with:

- valid:

  Clip data frame with validated, adjusted times.

- n_violated:

  Integer count of dropped clips.
