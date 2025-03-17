#!/bin/bash
java -jar /public/home/zhuzhenghang/picard.jar MergeVcfs \
          I=*.g.vcf.gz \
          O=LH20230111TAIL-seq.vcf.gz
