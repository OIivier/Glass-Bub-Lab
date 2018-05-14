plot(CI2, IBI2, '*')
hold
for q=1:119
    text(CI2(q),IBI2(q),num2str(q))
end
hold off