% Code by Yuxuan Xia

function in_gate = ellipsoidalGatingPoint(W,Point,model)

if isempty(W)
    in_gate = [];
    return;
end

d = 2;

in_gate = false(size(W,2),1);

H = model.measmodel.H(Point.xp);
%Take the extent into account when doing ellipsoidal gating
S = H*Point.Pp*H' + model.measmodel.R;
S = (S + S')/2;

nu = W - repmat(model.measmodel.h(Point.xp),[1,size(W,2)]);
dist= sum((inv(chol(S))'*nu).^2);

%Returns measurement indices inside the gate
in_gate(dist<model.gamma) = true;

end
