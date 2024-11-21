# Add a new analysis file -------------------------------------------------

wflow_open("analysis/read_data.Rmd")

# After opening a new analysis file, do the following:

# change: author: "Jens Daniel Müller"
# change: date:  "`r format(Sys.time(), '%d %B, %Y')`"

# include link to new html file in _site.yml

# Finally, rebuild, push and push again from R-Studio remaining files not taken care of by workflowr


# Repeated comments during work on the project ----------------------------

# to check impact of latest updates
wflow_build()

# commit regular changes (locally) and rebuild site
wflow_publish(here::here(
  "analysis",
  c(
    "index.Rmd",
    "read_data.Rmd",
    "temperature_indices.Rmd"
  )
),
message = "rebuild after code review")

# commit regular changes (locally) and rebuild site
wflow_publish(all = TRUE, message = "rebuild after code review")

# commit changes including _site.yml (locally) and rebuild site
wflow_publish(c("analysis/*Rmd"), message = "include G19 comparison", republish = TRUE)

# commit changes including _site.yml (locally) and rebuild site in the specified order
wflow_publish(here::here(
  "analysis",
  c(
    "index.Rmd",
    "read_data.Rmd",
    "enso_indices.Rmd"
  )
),
message = "rebuild after code update and new enso file",
republish = TRUE) #USE THIS ONE



# Push latest version to GitHub
wflow_git_push()
vgfroh
jens-daniel-mueller
