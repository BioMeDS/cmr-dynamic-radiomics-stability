# Run this script in the @MRXCAT_CMR_CINE subfolder of the MRXCAT directory
for snr in {5,10,20,30};
do
	perl -pe 's/^(SNR.*= )[\d.]+/${1}'$snr'/;' ../CINEpar.m >snr${snr}.m 
done
