using Plots
x=collect(1:10);
y=randn(10);
plotly()
plot(x,y,color="blue")
#and I have made some changes
gr()
plot(x,y,color="blue")
