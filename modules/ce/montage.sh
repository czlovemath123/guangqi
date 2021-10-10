python ce.py
f1='l'
f2='temp'
f3='kappa'
f4='v'
output='explosion'
postfix='.png'
for i in {1..40}
do
    l=$(printf $f1%05d$postfix $i)
    temp=$(printf $f2%05d$postfix $i)
    kappa=$(printf $f3%05d$postfix $i)
    v=$(printf $f4%05d$postfix $i)
    combine=$(printf $output%05d$postfix $i)
    montage $l $temp $kappa $v -tile 2x2 -geometry +2+2 $combine
done
eog $output'00001.png'
