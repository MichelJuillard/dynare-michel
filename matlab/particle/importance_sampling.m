function State_Particles = importance_sampling(StateMuPost,StateSqrtPPost,StateWeightsPost,numP)
[Xdim,Gsecond] = size(StateMuPost) ;  
u = rand(numP,1);
[Nc,comp] = histc(u, cumsum([0; StateWeightsPost]));    
State_Particles = zeros(Xdim,numP);
for k=1:Gsecond
  idx = bsxfun(@eq,comp,k*ones(size(comp))) ;
  State_Particles(:,idx) = StateSqrtPPost(:,:,k)*randn(Xdim,Nc(k));
end
State_Particles= State_Particles + StateMuPost(:,comp); 
    
