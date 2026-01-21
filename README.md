# Install on cluster

Copy spec and precalculated results to server

```{bash}
scp atlas_small-spec.RDS cahlmann@readinglab.cs.ucl.ac.uk:~/shinyTreelabel_data/.
scp atlas_small.duckdb cahlmann@readinglab.cs.ucl.ac.uk:~/shinyTreelabel_data/.
```

Update shinyTreelabel package

```{r}
remotes::install_github("const-ae/shinyTreelabel")
```

Lastly, restart the server running in the tmux session

```{bash}
Rscript ~/run_shiny_server.R
```
