# col=9 for npe
NF==9&&$1==11 {av += $col; a2 += $col^2; n++}
#NF==8&&$1==13 {av += $8; a2 += $8^2; n++}
END  {print av/n,"+/-", sqrt(((a2/n)-(av/n)^2)/n), n}
