# gfff2neo scripts

* `reload-database.sh` is run via a `cron` (once a month) to rebuild the graph database.
    - `0 0 1 * * ${HOME}/reload-database.sh`
* `backup.sh` is to be run via a `cron` to backup named volumes.
