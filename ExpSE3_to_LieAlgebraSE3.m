function d_err_d_rp_lie = ExpSE3_to_LieAlgebraSE3(d_err_d_rp, Tc0c1) 


  d_err_d_RPc0c1 = d_err_d_rp;
   RPc0c1 = Tc0c1;
 

  Rc0c1 = Tc0c1(1:3, 1:3);
  rotc0c1 = rodrigues(Rc0c1);

  d_err_d_rp_lie = d_err_d_RPc0c1 * Adj(RPc0c1) * rightJacobianSE3Decoupled(rotc0c1);
  



%   return d_err_d_rpv_lie;
  end