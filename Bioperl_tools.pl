#!/usr/bin/perl 
# use strict;
use Bio::DB::Fasta;
use Bio::Seq; 
use Bio::SeqIO;
use Spreadsheet::WriteExcel;
use LWP::UserAgent; 
use LWP::Simple;  

sub Get_the_PepID_From_Fasta {
    my $A = shift;

    my @List1;
    my $seqio_obj = Bio::SeqIO->new(-file => $A, 
                                 -format => 'fasta' );
    while($seq_obj = $seqio_obj->next_seq){
          push @List1, $seq_obj->display_id;
          push @List1, "\n";
      }

    my @name = split(/\./, $A);
    my $b = $name[0];
    print ("The File ($b\_ID.txt) has been built! \n");
    
    open(OUT, ">","$b\_ID.txt");
    print OUT @List1;
    close(OUT);
    
}
sub Count_the_Num_of_PepIDs_From_Fasta {
    my $A = shift;
    my @List1;
    my @ListPa;
    my @ListBg;
    my @ListZn;
    my @ListDm;
    my $Num_Pa = 0;
    my $Num_Bg = 0;
    my $Num_Zn = 0;
    my $Num_Dm = 0; 
    if ($A =~ m/(.*)_Seq/){
        $Name = $1;
    };  
    my $seqio_obj = Bio::SeqIO->new(-file => $A, 
                                 -format => 'fasta' );
    while($seq_obj = $seqio_obj->next_seq){
          push @List1, $seq_obj->display_id;
          push @List1, "\n";
          if ($seq_obj->display_id =~ m/Pa/){
              push @ListPa, $seq_obj->display_id;
              push @ListPa, "\n";
              $Num_Pa++;
          }elsif($seq_obj->display_id =~ m/BGER/){
              push @ListBg, $seq_obj->display_id;
              push @ListBg, "\n";
              $Num_Bg++;
          }elsif($seq_obj->display_id =~ m/Zn/){
              push @ListZn, $seq_obj->display_id;
              push @ListZn, "\n";
              $Num_Zn++;
          }elsif($seq_obj->display_id =~ m/Dm/){
              push @ListDm, $seq_obj->display_id;
              push @ListDm, "\n";
              $Num_Dm++;
          };                        
    }
    print "$Name\_Num_Pa = ",$Num_Pa,"\n";
    print "$Name\_Num_Bg = ",$Num_Bg,"\n";
    print "$Name\_Num_Zn = ",$Num_Zn,"\n";
    print "$Name\_Num_Dm = ",$Num_Dm,"\n";    
    push @List_Num, "Num_Pa = ",$Num_Pa,"\n";
    push @List_Num, "Num_Bg = ",$Num_Bg,"\n";
    push @List_Num, "Num_Zn = ",$Num_Zn,"\n";
    push @List_Num, "Num_Dm = ",$Num_Dm,"\n";
    
    open(OUT, ">$Name\_Num_Pa.txt");
    print OUT @ListPa;
    close(OUT);
    
    open(OUT, ">$Name\_Num_Bg.txt");
    print OUT @ListBg;
    close(OUT);
    
    open(OUT, ">$Name\_Num_Zn.txt");
    print OUT @ListZn;
    close(OUT);
    
    open(OUT, ">$Name\_Num_Dm.txt");
    print OUT @ListDm;
    close(OUT);
    
    open(OUT, ">$Name\_Num_Num.txt");
    print OUT @List_Num;
    close(OUT);
    
}
sub Sort_by_the_length_of_pep_2{
    my $c = shift;
    my $z = shift;
    my $seq;
    my %hash;
    my $n=1;
    my @List1;
    my @List2;
    
    # open(FD,"$c"); #
    # my @lines=<FD>;
    # my $m=0;
    # my $length = @lines;
    # close(FD);
    
    my $in=Bio::SeqIO->new(-file=>$c,-format=>"fasta");
    
    while($seq=$in->next_seq()){
        if (length($seq->seq) > $z){
            $hash{$seq->id}=length($seq->seq);
            push @List2, ">";
            push @List2, $seq->id;
            push @List2, "\n";
            push @List2, $seq->seq;
            push @List2, "\n";
            }
        }
    print @List2;
    # # foreach my $value ( sort values %hash ) {
    # foreach my $key (sort {$a cmp $b} keys %hash){       #Sort by key
    foreach my $key ( sort  { $hash{$a} <=> $hash{$b} } keys %hash ) {  #sort by value
         my $value = $hash{$key};
         print $n."\t".$key."\t".$value."\n";
         $n++;
         push @List1,"n";
         push @List1,"\t"; 
         push @List1, $key;
         push @List1,"\t"; 
         push @List1, $value;
         push @List1,"\n"; 
    };
    
    my @name = split(/\./, $c);
    my $d = $name[0];
    
    print ("The File ($d\_More_than_$z.txt) has been built! \n\n");
    print ("The File ($d\_More_than_$z.fasta) has been built! \n\n");
    
    open (OUT1,">","$d\_More_than_$z.txt");
    print OUT1 @List1;
    close(OUT1);
    
    open (OUT2,">","$d\_More_than_$z.fasta");
    print OUT2 @List2;
    close(OUT2);
};

sub Sort_by_the_DisplayID_of_pep{
    my $c = shift;
    
    my $seq;
    my %hash;
    my $n=1;
    my @List1;
    my $in=Bio::SeqIO->new(-file=>$c,
                           -format=>"fasta");
    while ($seq=$in->next_seq())
    {
       $hash{$seq->id}=$seq->seq(); #计算的是序列长度，序列的长度被存入hash表中
       # print $seq->id."\t".$seq->seq."\n";#$seq->seq()."\n";# 直接输入，输出的结果与上面awk的方法是一致的
    };
    
    # foreach my $value ( sort values %hash ) {
    foreach my $key (sort {$a cmp $b} keys %hash){       #Sort by key
    # foreach my $key ( sort  { $hash{$a} <=> $hash{$b} } keys %hash ) {  #sort by value
         my $value = $hash{$key};
         print $n."\t".$key."\t"."\n";
         # print $n."\t".$key."\t".$value."\n";
         $n++;
         push @List1,">";
         push @List1, $key;
         push @List1,"\n"; 
         push @List1, $value;
         push @List1,"\n"; 
    };

    my @name = split(/\./, $c);
    my $d = $name[0];
    print ("The File ($d\_Sorted.fasta) has been built! \n");
    
    open (OUT,">","$d\_Sorted.fasta");
    print OUT @List1;
    close(OUT);

};
sub Url_auto{
    my $file0 = shift; #query_pep_id
    my $file1 = shift; #Dm_pep_id
    my $file2 = shift; #Final_pep_id
    
    open(FD1,$file0);
    my @lines=<FD1>;
    my $length1 = @lines;
    close(FD1);

    open (FD2,$file1);        #
    my @array=<FD2>;
    my $length2 = @array;
    for($i=0;$i<$length1;$i++){
        my $url="http://202.127.22.44/tools/cockroach_gene.cgi?id=$lines[$i]";  
        my $content= get $url;  
        die"Couldn't get $url"unless defined $content;
        for($j=0;$j<$length2;$j++){
            my $id= $array[$j];
            chomp $id;
            # print $id; 
            if($content=~m/$id/){  
                push @List, $lines[$i];
                chomp $lines[$i];
                print "$lines[$i] \t It's ture! \n";  
            # }else{  
               #print"$lines[$i] Not ture.\n";  
            }; 
        };
    };
    close (FD2);
    
    my @FinalList = grep {++$count{$_} < 2} @List;
    print @FinalList;
    
    open (OUT1,">","$file2");   #
    print OUT1 @FinalList;
    close(OUT1);
}
sub RemoveDuplicates{
    my $c = shift;
    open File,"<",$c;
    my @List1 = <File>;
    my %count1;
    close File;
    my @List2 = grep {++$count{$_} < 2} @List1;
    my @List3 = sort { $a cmp $b } @List2;
    my @List4 = grep {$count{$_} > 1} @List1;
    my @List5 = grep {++$count1{$_} < 2} @List4;
    if (defined $List4[0]){
        print "Repetitive nameid Existes in $c!\n",@List5,"\n";
    }else{
        print ("All Nameids are Unique!\n");
    }
    my @path = split(/\./, $c);
    my $d = $path[0];
    print "The File ($d\_.Sort) has been built!\n";
    
    open  Out,">","$d\_.Sort";
    print Out @List3;
    close Out;
    open  Out,">","$d\_.repeat";
    print Out @List5;
    close Out;
    }
sub Get_Tags_from_NCBI_fasta{
    my $a = shift;
    my $b = shift;
    my $db = Bio::DB::Fasta->new($a);
    open FD,"<","$a";
    my @lines=<FD>;
    close FD;
    my $length=@lines;
    my @nameid;
    my @List;
    my @Finalfasta;
    $seqio_obj = Bio::SeqIO->new(-file => $a, 
                                 -format => 'fasta' );
    while($seq_obj = $seqio_obj->next_seq){
          push @List, $seq_obj->display_id.$seq_obj->desc;                                    
    }
    my $m = 0;
    foreach (@List){
          my $seqstring = $db->seq($_);
          push @Finalfasta,">Zn".$nameid[$m]."\n".$seqstring."\n";
          $m++;                     
    }
    open OUT,">","$b";
    print OUT @Finalfasta;
    close OUT;
};
sub Get_Tags_from_flybase_fasta{
    my $a = shift;
    my $db = Bio::DB::Fasta->new($a);
    open FD,"<","$a";
    my @lines=<FD>;
    close FD;
    my $length=@lines;
    my @nameid;
    my @List;
    my @Finalfasta;
    for($i=0;$i<$length;$i++){
          $name=$lines[$i];
          #if($name=~ m/Glutathione S-transferase (.*)/){ #Set the keywords
          # if($name=~ m/P450 (.*)/){ #Set the keywords
          # if($name=~ m/name=(.*); dbxref/){
              # if($name=~ m/name=(.*); parent=/){
           # if($name=~ m/ation_IDs:(\w+),/){
               # if($name=~ m/name=(.*); dbxref/){
             if($name=~ m/parent=(\w+),/){
                # print $1,"\n";
                push @nameid, $1;
                # push @nameid, "\n";
                };
    };
    # print @nameid, "\n";
    $seqio_obj = Bio::SeqIO->new(-file => $a, 
                                 -format => 'fasta' );
    while($seq_obj = $seqio_obj->next_seq){
          push @List, $seq_obj->display_id . "\n";             
          # print $seq_obj ->id,"\n";                          
    }
    my $m=0;
    foreach (@List){
          my $seqstring = $db->seq($_);
          push @Finalfasta,">Dm".$nameid[$m]."\n".$seqstring."\n";
          $m++;                      
    };
    my @name = split(/\./, $a);
    my $d = $name[0];
    print ("The File ($d\_Parent.fasta) has been built!\n");
    open OUT,">","$d\_Parent.fasta";
    print OUT @Finalfasta;
    close OUT;
    
    # open (OUT,">$b");
    # print OUT @List;
    # close(OUT);
    
    
};
sub Put_the_Pepid_into_Excel{
    my $a = shift;
    my $b = shift;    
    my $xl = Spreadsheet::WriteExcel->new("$b");
    my $xlsheet = $xl->add_worksheet("sheet1");  
    my $rptheader = $xl->add_format(); # Add a format
    $rptheader->set_bold();
    $rptheader->set_size('12');
    $rptheader->set_font('Century Gothic');
    my $normcell = $xl->add_format(); # Add a format
    $normcell->set_size('11');
    
    open(FD,"$a"); #
    my @lines=<FD>;
    my @List;
    my $m=0;
    my $length = @lines;
    close(FD);

    while($m<$length){
        my $id2=$lines[$m];
        # my $seq=$lines[$m+1];
        $m=$m+2;
        if ($id2=~m/>(.*)/){
            push @List, $1;
        };
    };
    my $num = 1;
    $length++;
    while ($num<$length){
        my @col_name = ("A".."ZZ");
        my $list_1 = shift @List;  
        my $col = shift @col_name;
        $xlsheet->write ("$col"."$num",$list_1,$normcell);
        $num++;
    };
}; 
sub Sort_by_the_length_of_pep{
    my $c = shift;
    my $seq;
    my %hash;
    my $n=1;
    my @List1;
    my $in=Bio::SeqIO->new(-file=>$c,-format=>"fasta");
    while ($seq=$in->next_seq()){
       $hash{$seq->id}=length($seq->seq()); #计算的是序列长度，序列的长度被存入hash表中
    }
    # foreach my $value ( sort values %hash ) {
    # foreach my $key (sort {$a cmp $b} keys %hash){       #Sort by key
    foreach my $key ( sort  { $hash{$a} <=> $hash{$b} } keys %hash ) {  #sort by value
         my $value = $hash{$key};
         print $n."\t".$key."\t".$value."\n";
         $n++;
         push @List1,$key."\t".$value."\n";
    }
    my @name = split(/\./, $c);
    my $d = $name[0];
    print ("The File ($d\_Sort_by_Length.txt) has been built! \n");
    open  OUT,">","$d\_Sort_by_Length.txt";
    print OUT @List1;
    close OUT;

};
sub Jiont_the_PepSeq_Files{
    my @PepList;
    my $length = @_;
    for($i=1;$i<$length;$i++){
        my $a = shift;
        $seqio_obj = Bio::SeqIO->new(-file => $a, 
                                     -format => 'fasta' );
        while($seq_obj = $seqio_obj->next_seq){
            push @PepList, ">".$seq_obj->display_id." ".$seq_obj->desc."\n".$seq_obj->seq."\n";              
            }
        };
        my $b = shift;
        open (OUT,">$b");
        print OUT @PepList;
        close(OUT);
    };
sub Get_Pepseq_by_ID{
    my $a = shift;
    my $b = shift;
    my $c = shift;
    my $db = Bio::DB::Fasta->new($a);
    open(FD,"$b"); 
    my @lines=<FD>;
    my @List;
    close(FD);

    foreach (@lines){
          chomp $_;
          my $seqstring = $db->seq($_);
          push @List,">";
          push @List, $_;
          push @List,"\n"; 
          push @List, $seqstring;
          push @List,"\n";                      
      };
    open (OUT,">$c");
    print OUT @List;
    close(OUT);
};

sub reverse_complement{
        my $dna = shift;

        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}                  

sub Split_the_PepSeq_Files{
    my $a = shift;
    my $b = shift;##path
    
    $seqio_obj = Bio::SeqIO->new(-file => $a, 
                                 -format => 'fasta' );
    while($seq_obj = $seqio_obj->next_seq){
      $id = $seq_obj->display_id;
      open OUT,">","$b/$id\_query\.fasta";
      open OUT1,">","$b/$id\_query\.id";
      print OUT ">".$seq_obj->display_id."\t".$seq_obj->desc."\n".$seq_obj->seq."\n";
      print OUT1 $seq_obj->display_id."\n";
      close OUT; 
      close OUT1;
      #push @PepList, ">".$seq_obj->display_id." ".$seq_obj->desc."\n".$seq_obj->seq."\n";              
      }
    }

Get_the_PepID_From_Fasta("");

# [Datebase_path] 
# $Datebase = "path";

# Get_Tags_from_flybase_fasta("$Input_name");
# Count_the_Num_of_PepIDs_From_Fasta("$Input_name");
# Sort_by_the_length_of_pep_2("$Input_name",200); 
# Sort_by_the_DisplayID_of_pep("$Input_name");
# Sort_by_the_length_of_pep("$Input_name");
# Jiont_the_PepSeq_Files("$Input_name_1",
#                        "$Input_name_2",
#                        "$Input_name_3",
#                        "$Input_name_4",
#                        "$Input_name_5",
#                        "$Output_name")

# Get_Pepseq_by_ID("$Datebase","$Input_name","$Output_name")
# Get_the_PepID_From_Fasta("$Input_name")
# RemoveDuplicates("$Input_name")
# Split_the_PepSeq_Files ("$Input_name","$path")
