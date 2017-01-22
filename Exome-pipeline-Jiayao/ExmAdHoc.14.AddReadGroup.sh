



java -jar picard.jar AddOrReplaceReadGroups \
      I=input.bam \
	        O=output.bam \
			      RGID=4 \
				        RGLB=lib1 \
						      RGPL=illumina \
							        RGPU=unit1 \
					RGSM=20
