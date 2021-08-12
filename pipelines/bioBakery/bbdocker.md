# Biobakery Workflows Docker

  * ec-2 configuration
    * can install docker there 
    * used Debian ami-628ad918
  * Instance that I created: biobakery (bcbb-aws)
    * add EBS volume to get mount point for additional storage
    * mount EBS volume <https://devopscube.com/mount-ebs-volume-ec2-instance/>
      * need to remount every time you start the instance
    * Set up AWS instance


```Bash
sudo apt-get install htop
sudo apt-get install emacs25
```

  * install docker
    * <https://docs.docker.com/engine/installation/linux/docker-ce/debian/>
    * post installation to run as non-sudo <https://docs.docker.com/engine/installation/linux/linux-postinstall/>



  * Change where docker stores images and volumes to the mounted EBS volume
    * stop the docker service before editing daemon.json, then restart
    ```bash
    service docker {start|stop|restart|status}
    ```
    * edit/create /etc/docker/daemon.json - replace graph path with EBS volume mountpoint
    ```json
    {
        "graph": "/mnt/docker",
        "storage-driver": "overlay"
    }
    ```

  * [biobakery workflows repository](https://hub.docker.com/r/biobakery/workflows/)
  * modify docker image permissions, etc.


```bash
docker pull biobakery/workflows
```

```Bash
# give sudo to default user linuxbrew
docker run -it --user root biobakery/workflows
apt-get install sudo
usermod -aG sudo linuxbrew
passwd linuxbrew #password is "nephele2"
exit

# commit changes to new dev branch
docker container ls -a #note container id
docker commit container-id biobakery/workflows:dev

# should have sudo now
docker run -it --user linuxbrew biobakery/workflows:dev

# install helpful stuff
sudo apt-get update
sudo apt-get install wget
sudo apt-get install htop
sudo apt-get install screen
sudo apt-get install python3-pip emacs
sudo pip3 install requests boto3 sh html2text sqlservice jinja2 humanize numpy biom-format
```


  * Running docker container after pulling <https://github.com/docker/labs/blob/master/beginner/chapters/alpine.md>
  * Modify a docker image <https://gist.github.com/glamp/74188691c91d52770807>


[Biobakery Workflow](https://bitbucket.org/biobakery/biobakery_workflows/wiki/Home)

  * download databases for WGS pipeline


```bash
# run docker interactively
docker run -it --user linuxbrew biobakery/workflows:dev

# mount /opt as volume so databases don't get committed - reduces image size
docker run -it -v /opt --user linuxbrew biobakery/workflows:dev

sudo /home/linuxbrew/.linuxbrew/bin/biobakery_workflows_databases --install wmgx --location <dir_name>
# these databases are 32G - installs to /opt/biobakery_workflows_databases

```


  * Bug: Human genome databases are in subfolder and wmgx kneaddata expects them to be in top folder
    * [Opened issue](https://bitbucket.org/biobakery/biobakery/issues/38/path-to-human-bowtie-index-is-incorrect-in)
    * Fix:


```bash
cd /opt/biobakery_workflows_databases/kneaddata_db_human_genome/
sudo ln -s Homo_sapiens_Bowtie2_v0.1/* .
```


  * Should commit changes to container after downloading dbs to save config changes made by biobakery_workflows_databases.  In case we forget to do that, can set config by the following:

    ```bash
    # if you forget to commit changes to container after downloading dbs, can change config:
    humann2_config --update database_folders utility_mapping /opt/biobakery_workflows_databases/humann2/utility_mapping/
    humann2_config --update database_folders protein /opt/biobakery_workflows_databases/humann2/uniref/
    humann2_config --update database_folders nucleotide /opt/biobakery_workflows_databases/humann2/chocophlan/
    ```

    ​

  * scp/rsync nephele2 git repo to /home/admin/git/nephele2 in instance

  * Run docker with git repo, database and data mounted


```bash
docker volume ls #note volume-id - this is where database is stored

docker run -it -v volume-id:/opt -v /home/admin/data:/mnt/EFS/user_uploads -v /home/admin/git/nephele2:/home/linuxbrew/nephele2 -v /home/admin/scripts:/home/linuxbrew/scripts --user linuxbrew biobakery/workflows:dev

```

* set PYTHONPATH in docker container .bashrc:

```bash
  export PYTHONPATH=$PYTHONPATH:/home/linuxbrew
```

* Appendix - Other cool stuff

  - Can provision docker image to AWS directly from laptop with [docker-machine](https://docs.docker.com/machine/examples/aws/)

* Useful commands

  - docker top - on running container
  - docker info - shows number of running containers, images, etc.
  - docker start container-id - starts stopped container
  - docker attach container-id - attaches running container (have to do Ctrl-C to get bash prompt? )
  - docker ps - lists running containers; docker ps -a - lists all containers, even those that are stopped
  - docker system df - lists docker disk usage
  - docker kill container-id - kills running container -- don't forget to do docker system df to see if the stopped container is using up space - may have to do `docker prune container-id`
  - [Leave container bash prompt without stopping container](https://stackoverflow.com/a/36791612) - Ctrl+p+q  
  - docker stats - collect cpu usage, file i/o, memory info; [Simple profiling script](https://github.niaid.nih.gov/subramanianp4/dotfiles/blob/master/aws/dockerstats.sh)


* Update biobakery/workflows docker container image from 0.3.1 to 0.13.2
  - Create a new AMI (with Packer or simply take the EC2 snapshot) that has Docker image with this Dockerfile
  ```bash
  FROM biobakery/workflows:0.13.2

  RUN humann2_config --update database_folders utility_mapping /opt/biobakery_workflows_databases/humann2/utility_mapping/ &&\
      humann2_config --update database_folders protein /opt/biobakery_workflows_databases/humann2/uniref/ &&\
      humann2_config --update database_folders nucleotide /opt/biobakery_workflows_databases/humann2/chocophlan/
  ```
  - Share AMI permission between DEV, QA and PROD
  - Download database for biobakery/workflows
    - Create a new folder for biobakery/workflows's database /mnt/EFS/dbs/biobakery_workflows_databases_0.13.2
    - Mount EFS to the machine that you are running on
    - Run the command
    ```bash
    sudo docker run --mount type=bind,source=/mnt/EFS/dbs/biobakery_workflows_databases_0.13.2,target=/mnt/EFS/dbs/biobakery_workflows_databases_0.13.2 --rm biobakery/workflows:0.13.2 biobakery_workflows_databases --install wmgx --location /mnt/EFS/dbs/biobakery_workflows_databases_0.13.2
    ```
