M$Restrictions = {
   DU[i_, j_] :> 0 /; NumericQ[i] && NumericQ[j] && (i!= j),
   DU[i_, i_] :> yuku[i],
   DD[i_, j_] :> 0 /; NumericQ[i] && NumericQ[j] && (i!= j),
   DD[i_, i_] :> yukd[i],
   DL[i_, j_] :> 0 /; NumericQ[i] && NumericQ[j] && (i!= j),
   DL[i_, i_] :> yukl[i],
   CKM[i_Integer, 3] :> 0 /; i != 3,
   CKM[3, i_Integer] :> 0 /; i != 3,
   CKM[3, 3] :> 1
}
