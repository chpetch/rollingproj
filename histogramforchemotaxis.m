ccolor = {'b','b','b','r','r','r','r','g','g','g','g'};
numgauss = [2,2,2,2,2,2,2,1,1,1,1];
A = {HW1_1,HW2_1,XW_1, ...
    HWTNFa1_1,HWTNFa2_1,XWTNFa1_1,XWTNFa2_1, ...
    HW1_0,HW2_0,XW_0,HWTNFa1_0};
fname = {'mix'};
name = {'HW1_1_120min','HW2_1_120min','XW_1_120min', ...
    'HWTNFa1_1_120min','HWTNFa2_1_120min','XWTNFa1_1_120min','XWTNFa2_1_120min',...
    'HW1_0min','HW2_0min','XW_0min','HWTNFa1_0min'};
fldname = 'C:\Users\Chayakorn\OneDrive for Business\drhou\collected results\chemotaxis\17082015';
clf
for i = 1:1:11
    A{i}(find(A{i}==-1)) = [];
    [histFreq, histXout] = hist(A{i},[0:10:1200]);
    gcurve = gmdistribution.fit(A{i},numgauss(i));
    for j = 1:1:4
        figure(j)
        if j == 1
            hold on
            bar(histXout, histFreq/sum(histFreq)*100,ccolor{i})
            hold on
            plot(histXout, pdf(gcurve, histXout')*100,[ccolor{i} '-']);
        elseif j == 2
            plot(histXout, pdf(gcurve, histXout')*100,[ccolor{i} '-']);
            title(name{i})
            saveas(figure(2),[fldname,'\',name{i},'_.jpg'])
        elseif j == 3
            if i<8
                plot(histXout, pdf(gcurve, histXout')*100,[ccolor{i} '-']);
                hold on
                title('All curve')
            end
        elseif j == 4
            [histFreq, histXout] = hist(A{i},[0:10:1200]);
            hold on
            histdata = histFreq/sum(histFreq)*100;
            bar(histXout, histdata,ccolor{i})
            hold on
            pdfdata = pdf(gcurve, histXout')*100;
            
            mpdfdata = pdfdata/max(pdfdata);
            plot(histXout, mpdfdata*max(histdata),[ccolor{i} '-']);
        end
    end
    result.pdf{i} = pdfdata;
    result.bin{i} = histXout;
    result.mean{i} = gcurve.mu;
    result.max{i} = max(histdata);
end
saveas(figure(1),[fldname,'\',fname{1},'_histandcurve.jpg'])
saveas(figure(3),[fldname,'\',fname{1},'_curve0-120.jpg'])
saveas(figure(4),[fldname,'\',fname{1},'_histandscaledcurve.jpg'])