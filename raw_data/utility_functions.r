# automatically convert text to something easy to work with
ensnakeify <- function(x) {
  x %>%
    iconv(to="ASCII//TRANSLIT") %>% # remove accents
    str_replace_na() %>% # convert NA to string
    str_to_lower() %>% # convert to lower case
    str_replace_all(pattern="[^[:alnum:]]", replacement=" ") %>% # convert non-alphanumeric to space
    str_trim() %>% # trim leading and trailing spaces
    str_replace_all(pattern="\\s+", replacement="_") # convert remaining spaces to underscore
}

#
autosnake <- function(df){ # to use in pipe
  names(df) <- names(df) %>%
    iconv(to="ASCII//TRANSLIT") %>% # remove accents
    str_replace_na() %>% # convert NA to string
    str_to_lower() %>% # convert to lower case
    str_replace_all(pattern="[^[:alnum:]]", replacement=" ") %>% # convert non-alphanumeric to space
    str_trim() %>% # trim leading and trailing spaces
    str_replace_all(pattern="\\s+", replacement="_") # convert remaining spaces to underscore
  df
}

# rounding function that works on dates
round_any = function(x, accuracy, f=round){
  if (is.Date(x)){
    switch(f,
           round=round_date(x, accuracy),
           ceiling=ceiling_date(x, accuracy),
           floor=floor_date(x, accuracy)
    )
  } else {
    f( x / accuracy) * accuracy
  }
}
