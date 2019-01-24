FROM debian:stretch

MAINTAINER Dr. Asis Hallab <asis.hallab@gmail.com>

# update debian package list
RUN apt-get update

# install R and java
RUN apt-get install -y r-base openjdk-8-jre openjdk-8-jdk

# install personal packages
RUN apt-get install -y ssh neovim zsh tmux python-pip virtualenvwrapper \
	curl sudo locales locate less wget git 

# install indispensable packages
RUN apt-get install -y mcl fasttree mafft ncbi-blast+ wget hyphy-mpi \
	libssl-dev libcurl4-openssl-dev libxml2-dev exuberant-ctags

# add Debian sid repo 
RUN echo 'deb http://deb.debian.org/debian sid main' >>  /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y diamond-aligner/unstable paml/unstable
RUN apt-get clean


# install CAFE
RUN wget -O /opt/cafe https://downloads.sourceforge.net/project/cafehahnlab/cafe.linux.x86_64  
RUN chmod +x /opt/cafe 
RUN cp /opt/cafe /usr/local/bin/

# download, compile and install Blat
RUN wget -O /opt/blatSrc35.zip https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
RUN unzip /opt/blatSrc35.zip -d /opt/blatSrc

ENV MACHTYPE x86_64
RUN mkdir -p /root/bin/$MACHTYPE/ 
RUN make -C /opt/blatSrc/blatSrc
RUN cp /root/bin/$MACHTYPE/* /usr/local/bin/
#RUN rm -rf /root/bin

# install and configure locales
RUN apt-get install -y locales
RUN sed -i 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    echo 'LANG="en_US.UTF-8"'>/etc/default/locale

RUN dpkg-reconfigure --frontend=noninteractive locales
RUN update-locale LANG=en_US.UTF-8

# install R packages
RUN sh -c 'wget -O /home/tlaloc/install_GeneFamilies.R https://raw.githubusercontent.com/asishallab/GeneFamilies/master/exec/install_script.R'
RUN sh -c 'R -f /home/tlaloc/install_GeneFamilies.R'
RUN echo 'install.packages("roxygen2", repos="http://cran.us.r-project.org"); install.packages("formatR", repos="http://cran.us.r-project.org");' > /home/tlaloc/install_R_pkgs.R
RUN sh -c 'R -f /home/tlaloc/install_R_pkgs.R'
RUN sh -c 'rm /home/tlaloc/install_GeneFamilies.R'
