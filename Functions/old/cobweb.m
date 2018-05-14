function h=cobweb(mapfn,a,b,xo,niter)
% Draws a cobweb diagram of function mapfn, with domain [a,b]
% and starting value xo, with niter iterations.
% Returns handle to the diagram.
%
% The dynamical mapping function is plotted with a solid, blue line.
% The line y=x is plotted with a dashed, black line.
% The initial leg of the cobweb is plotted with a red dotted line.
% The main part of the cobweb is plotted with a red dot-dashed line,
% with arrows indicating the forward direction of the iteration.
%
% Input arguments:
%		mapfn is the function to be cobwebbed.
%			It can either be an explicit function of x,
%			or a MATLAB function M-file.
%		a is the lower limit of the domain
%		b is the upper limit of the domain
% 		xo is the starting value of x
%		niter is the number of iterations to compute and plot
%
% Examples:
%	To generate a cobweb diagram for a two worker TSS line with
%		r1=ratio of speed of first worker to second worker = 0.8,
%		and perturbed handoff position x2=0.35,
%		you can use either:
%
%			cobweb('tss2map',0,1,0.35,8)
%			cobweb('min(0.8*(1-x),1)',0,1,0.35,8)
%
%		where tss2map is a MATLAB function also available at the
%		course web page.  (r1 is set within that function M-file.)
%
%	To generate a cobweb diagram for a two worker TSS line with
%		r1=ratio of speed of first worker to second worker = 1.2,
%		and perturbed handoff position x2=0.5,
%		you can use either:
%
%			cobweb('tss2map',0,1,0.5,16)
%			cobweb('min(1.2*(1-x),1)',0,1,0.5,16)
%
%		where tss2map is a MATLAB function also available at the
%		course web page.  (r1 is set within that function M-file.)
%
% The number of iterations (last argument, niter), is best set by
%	experimenting until you get a nice-looking cobweb.


hold off
fplot(fcnchk(mapfn),[a b],[],[],'b-') % plot the dynamical mapping function
axis equal
axis([a b a b])
hold on
fplot('1*x',[a b],[],[],'k--') % plot the line y=x for reference
x=zeros(2*niter);
y=zeros(2*niter);
x(1)=xo;
y(1)=a;
plot(x(1),y(1),'ro')
x(2)=xo;
y(2)=feval(fcnchk(mapfn),xo);
if (y(2) > b) | (y(2) < a)
   error('Domain of map incorrectly specified.')
end
yarrow=(a+y(2))/2;
plot(xo,yarrow,'r^') % plot an up arrow
lrarrow(1)=0;
udarrow(1)=1;
for i=2:niter
   xold=x(2*(i-1));
   yold=y(2*(i-1));
   ynew=feval(fcnchk(mapfn),yold);
   x(2*i-1)=yold;
   y(2*i-1)=yold;
   lrarrow=sign(yold-xold);
   xarrow=(xold+yold)/2;
   if (lrarrow > 0) % plot right arrow
      plot(xarrow,yold,'r>')
   elseif (lrarrow < 0)
         plot(xarrow,yold,'r<') % plot left arrow
   end
   x(2*i)=yold;
   y(2*i)=ynew;
   if (ynew > b) | (ynew < a)
      error('Domain of map incorrectly specified.')
   end
   udarrow=sign(ynew-yold);
   yarrow=(ynew+yold)/2;
   if (udarrow > 0) % plot up arrow
      plot(yold,yarrow,'r^')
   elseif (udarrow < 0) % plot down arrow
         plot(yold,yarrow,'rV')
   end
end
plot(x(1:2),y(1:2),'r:') % plot initial leg with a dotted line
plot(x(2:2*niter),y(2:2*niter),'r-.')
	  % plot the rest of the iterations with dot-dashed line
plot(x(2*niter),y(2*niter),'rx')
xlabel('x')
ylabel('y')
title('Cobweb diagram')
h=figure(gcf);
ynew
return
