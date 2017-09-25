NF==1 {
    a_or=17.925; b_or=-0.057681;
    a_ex=13.753; b_ex=-0.046095;
    alor=exp(a_or+b_or*$1);
    alex=exp(a_ex+b_ex*$1);
    alfa=(alor+alex)/2.;
    print $1,1./alfa;
}
NF==3 {print $1,2/($2+$3)}
