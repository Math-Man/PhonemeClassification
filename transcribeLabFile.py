#Program is only going to create lab for this file name, in the segments
AudioFileName = "VOA_070206_0630";



segmentsReader = open("segments","r");
phonesReader = open("phones.txt","r");
labelFile = open("labelFile06BASIC.lab", "w+");
labelFilewSilence = open("labelFile06DETAILED.lab", "w+");
#line = segments.readline();
#words = line.split(" ");
#print(words[1]);

#Phoneme dictionary
phonemesList = [];
phoLine = phonesReader.readline();
while (phoLine):
    phoWords = phoLine.split("  ");
    phonemeWords = phoWords[0].rstrip().split("_");
    phonemesList.append({ 'phoneme':phoWords[0].rstrip(), 'index':phoWords[1].rstrip()});
    phoLine = phonesReader.readline();
print(phonemesList)
#print(phonemesList[0].get('phoneme'))
#a = next(item for item in phonemesList if item['index'] == '1')
#print(a)




segLine = segmentsReader.readline();
while (segLine):
    segWords = segLine.split(" ");
    if((segWords[1] is not None) and (segWords[1] == AudioFileName)): #Check if the file name matches
        #scan alignments
        alignmentsReader = open("alignments","r");
        alignLine = alignmentsReader.readline();
        while(alignLine):
            alignWords = alignLine.split(" ");
            if((alignWords[0] == segWords[0])): #if the speaker segments match

                #Take segment offset
                offset = segWords[2];
                #print(alignLine + " " + offset + " " + segWords[1]);
                for i in phonemesList:
                    ci = (float(i.get('index')))
                    cc = (float(alignWords[4]));
                    if(ci == cc):
                        currentPhoneme = (i.get('phoneme'));
                        #print(currentPhoneme)

                #print(currentPhoneme)
                #Write start+offset, end+offset, phoneme(only word)
                labelFilewSilence.write( str( '%.5f' % (float(offset) + float(alignWords[2])) ) + " " + str( '%.5f' % (float(offset) + float(alignWords[2]) + float(alignWords[3])) ) + " " + currentPhoneme + "\n");

                #remove silences and after _
                scrapedPhone = currentPhoneme.split("_");
                scrapedPhone = scrapedPhone[0];
                if(scrapedPhone != 'SPN'): # (scrapedPhone != 'SIL') &
                    labelFile.write( str( '%.5f' % (float(offset) + float(alignWords[2])) ) + " " + str( '%.5f' % (float(offset) + float(alignWords[2]) + float(alignWords[3])) ) + " " + scrapedPhone + "\n");


            alignLine = alignmentsReader.readline();
        alignmentsReader.close();
    segLine = segmentsReader.readline();
segmentsReader.close();

print("done")
