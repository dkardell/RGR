function [B,N,fn] = build_bnd(Param,Grid,I) 

% author: Dominik A. Kardell
% modified from: Evan J. Ramos and Marc A. Hesse
% date: December 2020

% Description: 
% This function computes the operators and r.h.s vectors for both Dirichlet 
% and Neumann boundary conditions. 

% Input: 
% Grid = structure containing all pertinent information about the grid. 
% Param = structure containing all information about the physical problem 
% in particular this function needs the fields 
% Param.dof_dir = Nc by 1 column vector containing 
% the dofs of the Dirichlet boundary. 
% Param.dof_neu = N by 1 column vector containing 
% the dofs of the Neumann boundary. 
% Param.qb = prescribed fluxes on Neuman bnd.  

% Output: 
% B = Nc by N matrix of the Dirichlet constraints 
% N = (N-Nc) by (N-Nc) matrix of the nullspace of B 
% fn = N by 1 r.h.s. vector of Neuman contributions
    
    B = I(Param.dof_dir,:);
    N = I;
    N(:,Param.dof_dir) = [];
    
    if ~isempty(Param.qb)
        dx = Grid.V(Param.dof_neu)./Grid.A(Param.dof_f_neu);
        fn = spalloc(Grid.N,1,length(Param.dof_neu));
        fn(Param.dof_neu) = Param.qb./dx;
        
        if length(Param.dof_neu) > 1 && Grid.Ny == 1
            Param.dof_neu(end) = -Param.dof_neu(end);
        end
        
    else
        
        fn = sparse(Grid.N,1);
        
    end
    
end