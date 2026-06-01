code_dir <- "Codes"
if (!dir.exists(code_dir)) {
  stop("Codes/ directory not found. Run this script from the project root.")
}

r_files <- list.files(
  code_dir,
  pattern = "\\.R$",
  recursive = FALSE,
  full.names = TRUE
)

parse_errors <- list()
for (file in r_files) {
  ok <- tryCatch(
    {
      parse(file)
      TRUE
    },
    error = function(err) {
      parse_errors[[file]] <<- conditionMessage(err)
      FALSE
    }
  )

  message(if (ok) "OK  " else "BAD ", file)
}

if (length(parse_errors) > 0) {
  print(parse_errors)
  stop("One or more active R scripts failed to parse.")
}

message(length(r_files), " active R scripts parse OK.")
