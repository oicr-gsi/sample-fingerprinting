#bin/bash
set -e -o pipefail

repository="git@github.com:oicr-gsi/sample-fingerprinting.git"
localFolder="/home/ubuntu/sample-fingerprinting"

if cd $localFolder; then
  git pull
else 
  git clone "$repository" "$localFolder"
fi

git checkout TNG_pipeline_day

docker run -it --rm \
  -v "${PWD}":"${localFolder}" \
  -v "${HOME}/.m2":"/root/.m2" \
  -w "${localFolder}" \
  maven:3.3-jdk-7 \
  mvn clean install

cd workflow-fingerprint-collector/target

mkdir -p workflows && mkdir -p datastore && mkdir -p packaged-bundles
chmod a+wrx workflows && chmod a+wrx datastore && chmod a+wrx packaged-bundles
mkdir -p /home/ubuntu/zipped_wf

# doesn't work due to broken wrench (not seqwaremaven) reference
# docker run --rm -h master -t -v `pwd`/datastore:/mnt/datastore -v "${PWD}":"/home/seqware" -v "${HOME}/packaged-bundles":"/home/seqware/packaged-bundles" -i seqware/seqware_whitestar:1.1.2-java8 seqware bundle package --dir /home/seqware/Workflow_Bundle_FingerprintCollector_1.2.1_SeqWare_1.1.0 --to /home/seqware/packaged-bundles --no-metadata


echo "successfully buil"
