function [Q,details] = fig3_modes(M,t,goTime)
cr=M(:,:,1); cl=M(:,:,2); er=M(:,:,3); el=M(:,:,4);
avg=(cr+cl)/2; A=zeros(size(cr,1),8);
A(:,1)=meanwindow((cr+er-cl-el)/2,t,0,1.3);
A(:,2)=meanwindow((cr+el-cl-er)/2,t,1.3,goTime);
A(:,3)=meanwindow((cr+cl-er-el)/2,t,goTime,goTime+1.3);
A(:,4)=meanwindow(cr-cl,t,.2,.4);
A(:,5)=meanwindow(cr-cl,t,goTime-.3,goTime-.1);
A(:,6)=meanwindow(cr-cl,t,goTime+.1,goTime+.3);
A(:,7)=meanwindow(avg,t,goTime-.3,goTime-.1)-meanwindow(avg,t,-.3,-.1);
A(:,8)=meanwindow(avg,t,goTime,goTime+.1)-meanwindow(avg,t,goTime-.1,goTime);
lengths=sqrt(sum(A.^2,1));
assert(all(lengths>eps),'A task-defined mode is zero');
A=A./lengths;
assert(rank(A)==8,'Task-defined modes are dependent; original mode indices would be ambiguous');
Q=zeros(size(A));
for i=1:8
    u=A(:,i);v=zeros(size(u));
    for j=1:i-1,v=v-dot(u,Q(:,j))*Q(:,j);end
    w=u+v;Q(:,i)=w/norm(w);
end
assert(norm(Q'*Q-eye(8),'fro')<1e-8,'Mode orthogonality check failed');
details.rawNormalizedModes=A;
details.order={'stimulus','choice','outcome','early_sample','late_delay','early_response','ramp','go'};
details.fit='All available trials, descriptive in-sample modes; no held-out validation';
end
function y=meanwindow(X,t,a,b)
sel=t>a & t<b;assert(any(sel),'Empty mode window');y=mean(X(:,sel),2);
end
