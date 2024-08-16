import csv
import statistics as statistics

#write to file
with open('latex_table.txt', 'w') as f:
    f.write("\hline \n")
    #f.writelines(["Method & Min. & Med. & Avg. & Max. & Std. & Win & Fail \\\ ", "\n"])
    f.writelines(["Method & Med. & Avg. & Std. & Win & Fail \\\ ", "\n"])
    f.write("\hline \n")

    for K in [10,100,500,1000,5000,10000]:
        
        #opening the CSV file
        with open(str(K)+'.csv', mode ='r')as file:
            
            #reading the CSV file
            csvFile = list(csv.reader(file))
            methods = csvFile[0]
            
            #List of statistics
            min_list = []
            max_list = []
            med_list = []
            avg_list = []
            std_list = []
            win_list = []
            fail_list = []
            
            #Get statistics
            for column in range(len(methods)):
                results = [float(csvFile[i][column]) for i in range(1,len(csvFile))]
                max_list.append(max(results))
                min_list.append(min(results))
                avg_list.append(statistics.mean(results))
                med_list.append(statistics.median(results))
                std_list.append(statistics.stdev(results))
                fail_list.append(sum(x > 1800 for x in results))
                win_list.append(sum(
                    [ all( [ float(csvFile[i][j]) >= float(csvFile[i][column]) 
                            and float(csvFile[i][column]) < 1800
                            for j in range(len(methods)) ] )
                     for i in range(1,len(csvFile))]))
            
    
            #Write heading
            dic = {10:"a) 10", 100:"b) 100", 500:"c) 500", 1000:"d) 1,000", 5000:"e) 5,000", 10000:"f) 10,000" }
            #f.writelines("\multicolumn{8}{c}{" + dic[K] + " segments} \\\ \n")
            f.writelines("\multicolumn{6}{c}{" + dic[K] + " segments} \\\ \n")
            
            
            #Sort methods after median
            sort_med = med_list.copy()
            sort_med.sort()
            
            #Method dictionary
            method_dic = {"sBB" : "sBB", "DisaggLogarithmic" : "DLog", "ZigZagInteger" : "ZZI", 
                          "Logarithmic" : "Log", "ZigZag" : "ZZB"}
            
            #Writing results to .txt
            if K in [10,100]: formatstring = '{0:,.2f}'
            elif K in [500,1000,5000]: formatstring = '{0:,.1f}'
            elif K in [10000]: formatstring = '{0:,.0f}'
            
            for val in sort_med:
                i = med_list.index(val)
                f.writelines(method_dic[methods[i]]+ " & " +
                             #'{0:,.2f}'.format((round(min_list[i],2)))+ " & " +
                             "\\textbf{" + formatstring.format(round(med_list[i],2))+ "}" + " & " +
                             formatstring.format(round(avg_list[i],2))+ " & " +
                             #'{0:,.2f}'.format((round(max_list[i],2)))+ " & " +
                             formatstring.format(round(std_list[i],2))+ " & " +
                             str(win_list[i])+ " & " +
                             str(fail_list[i])+
                             "\\\ \n")
            
    f.write("\hline \n")

        
        
