function a=test(eps)
    a=zeros(6,2);
    a(1,1)=eval_interpolator_c(1, eps);
    a(2,1)=eval_interpolator_c(2, eps);
    a(3,1)=eval_interpolator_c(3, eps);
    a(4,1)=eval_interpolator_c(4, eps);
    a(5,1)=eval_interpolator_c(5, eps);
    a(6,1)=eval_interpolator_c(6, eps);
    
    a(1,2)=eval_interpolator_d(1, eps);
    a(2,2)=eval_interpolator_d(2, eps);
    a(3,2)=eval_interpolator_d(3, eps);
    a(4,2)=eval_interpolator_d(4, eps);
    a(5,2)=eval_interpolator_d(5, eps);
    a(6,2)=eval_interpolator_d(6, eps);
end

