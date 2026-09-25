# syntax=docker/dockerfile:1
FROM dfam/tetools:latest

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN <<'EOF'
set -eux
export DEBIAN_FRONTEND=noninteractive

apt-get -y update
apt-get install --assume-yes software-properties-common
# add-apt-repository universe
apt-get update
apt-get install --assume-yes python3-pip git build-essential zlib1g-dev libcereal-dev libjellyfish-2.0-dev pkg-config cmake r-base-core gawk autoconf pigz rustc cargo

apt-get -y install \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    curl libgomp1 \
    perl \
    libfile-which-perl \
    libtext-soundex-perl \
    libjson-perl liburi-perl libwww-perl \
    libdevel-size-perl \
    bedtools \
    ncbi-blast+

apt-get install --assume-yes tabix libbz2-dev liblzma-dev libgsl-dev libperl-dev bzip2
rm -rf /var/lib/apt/lists/*
EOF


RUN <<'EOF'
set -eux
cd "${HOME}"
git clone https://github.com/marbl/Winnowmap.git
cd Winnowmap
make -j8
cp bin/* /usr/local/bin/
cd ..
rm -r Winnowmap
EOF

RUN <<'EOF'
set -eux
cd "${HOME}"
git clone --recursive https://github.com/samtools/htslib.git
cd htslib
autoreconf -i
./configure
make
make install
cd "${HOME}"

git clone https://github.com/samtools/samtools.git
cd samtools
autoheader
autoconf -Wno-syntax
./configure --without-curses
make
make install
cd "${HOME}"
rm -rf samtools

git clone https://github.com/samtools/bcftools.git
cd bcftools
autoheader
autoconf
./configure --enable-libgsl --enable-perl-filters
make
make install

cd "${HOME}"
rm -rf htslib
rm -rf bcftools
EOF

RUN <<'EOF'
set -eux
cd "${HOME}"
git clone https://github.com/fritzsedlazeck/SURVIVOR.git
cd SURVIVOR/Debug
make
cp SURVIVOR /usr/local/bin
cd "${HOME}"
rm -rf SURVIVOR

git clone https://github.com/lh3/minimap2
cd minimap2
make
cp minimap2 /usr/local/bin
cd "${HOME}"
rm -rf minimap2

git clone https://github.com/lh3/minigraph
cd minigraph
make
cp minigraph /usr/local/bin
cd "${HOME}"
rm -rf minigraph

# ULTRA finds the tandem repeats that go into total_repeat_span
# (bin/repmask_vcf.sh calls `ultra`). Check the tag against the published image.
git clone --branch v1.0.0 --depth 1 https://github.com/TravisWheelerLab/ULTRA.git
cd ULTRA
cmake .
make
cp ultra /usr/local/bin
cd "${HOME}"
rm -rf ULTRA
EOF

RUN <<'EOF'
set -eux

cd /tmp
git clone --branch v1.3.2 --depth 1 https://github.com/USCiLab/cereal.git
cd cereal
mkdir build && cd build
cmake -DJUST_INSTALL_CEREAL=ON ..
make install
cd /tmp && rm -rf cereal

mkdir /metadata
dpkg -l | grep jellyfish | tr -s " " | cut -d " " -f 2,3 > /metadata/jellyfish.lib.version
mkdir /repos
cd /repos
git clone https://github.com/eblerjana/pangenie.git
cd pangenie
mkdir build
cd build
cmake ..
make -j 4
cp src/PanGenie /usr/local/bin
cp src/PanGenie-index /usr/local/bin
cd ..
git rev-parse --short HEAD > /metadata/pangenie.git.version
cd "${HOME}"
rm -rf /repos/pangenie
EOF

RUN <<'EOF'
set -eux
#pip3 install --break-system-packages pyabpoa
pip3 install --break-system-packages pysam pyparsing svim-asm pandas polars vcfpy sniffles cigar truvari pyfaidx h5py
EOF

RUN <<'EOF'
R --slave -e 'install.packages(c("XML", "dplyr", "stringr", "tidyr", "readr", "vcfR", "optparse"), repos="https://cloud.r-project.org/")'
EOF

RUN <<'EOF'
set -eux
export DEBIAN_FRONTEND=noninteractive
cd "${HOME}"
# Install dependencies and some basic utilities.
apt-get -y update
apt-get -y install \
    aptitude \
    libgomp1 \
    perl \
    libfile-which-perl \
    libtext-soundex-perl \
    libjson-perl liburi-perl libwww-perl \
    libdevel-size-perl
aptitude install -y ~pstandard ~prequired \
    curl wget \
    vim nano \
    procps strace \
    libpam-systemd-

echo "PS1='(dfam-tetools \$(pwd))\\\$ '" >> /etc/bash.bashrc

apt-get -y install bc
apt-get clean --assume-yes
rm -rf /var/lib/apt/lists/*
EOF

RUN <<'EOF'
set -eux
wget -O /usr/local/bin/vg https://github.com/vgteam/vg/releases/download/v1.77.0/vg
chmod +x /usr/local/bin/vg

# pypy3 at /opt/pypy3: bin/subset_gaf.py filters every alignment line of every
# sample and its shebang is /opt/pypy3/bin/pypy3.
wget -qO /tmp/pypy3.tar.bz2 https://downloads.python.org/pypy/pypy3.10-v7.3.17-linux64.tar.bz2
mkdir -p /opt/pypy3
tar -xj -f /tmp/pypy3.tar.bz2 --strip-components=1 -C /opt/pypy3
rm /tmp/pypy3.tar.bz2
EOF


RUN <<'EOF'
set -eux
wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O "${HOME}/miniforge.sh"
bash "${HOME}/miniforge.sh" -b -p "${HOME}/miniforge"
rm -f "${HOME}/miniforge.sh"
"${HOME}/miniforge/bin/conda" create -y -n ga -c conda-forge -c bioconda graphaligner
cp "${HOME}/miniforge/envs/ga/bin/GraphAligner" /usr/local/bin/
"${HOME}/miniforge/envs/ga/bin/GraphAligner" --version
rm -rf "${HOME}/miniforge"
EOF

RUN <<'EOF'
set -eux
curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain 1.85.0
export PATH="${HOME}/.cargo/bin:${PATH}"
cd "${HOME}"
git clone https://github.com/cgroza/panmethyl
cd panmethyl/tagtobed
cargo build --release
cp \
    target/release/lift_mods \
    target/release/lift_offsets \
    target/release/lift_edges \
    target/release/tagtobed \
    /usr/local/bin/
cd "${HOME}"
rm -rf panmethyl
EOF

RUN <<'EOF'
set -eux
export DEBIAN_FRONTEND=noninteractive
apt-get remove --assume-yes git software-properties-common cmake make pkg-config build-essential autoconf
apt-get autoremove --assume-yes
apt-get clean --assume-yes
rm -rf /var/lib/apt/lists/*

command -v minigraph
command -v lift_mods
command -v lift_offsets
command -v lift_edges
python3 -c 'import polars'
EOF

ENV LC_ALL=C \
    LANG=C \
    PYTHONIOENCODING=utf8 \
    PATH=/opt/RepeatMasker:/opt/RepeatMasker/util:/opt/RepeatModeler:/opt/RepeatModeler/util:/opt/coseg:/opt/ucsc_tools:/opt:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
