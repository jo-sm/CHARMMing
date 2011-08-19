#!/usr/bin/perl
use POSIX;
$num=$ARGV[0];
$upload_dir = $ARGV[1];
$filename=$ARGV[2];
$upload_dir=~s/"//g;
$upload_dir=~tr/A-Z/a-z/;
#$upload_dir.="\/";
$filename=~s/"//g;
$filename=~tr/A-Z/a-z/;
$num=ceil($num);

#contains prime numbers between 5 and 100
@primes = qw(7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127 131);

#checks if a number mod another number is equal to 0
sub isprime
{
 $_[0] % $_[1];
}


$i=0;
$primeindex=@primes;
while(($i+1)<$primeindex && ($num > $primes[$i])) #while i is less than the size of the @primes array and while the $num is less than the prime num
{
 $value = $num % 2;
 if(&isprime($num,$primes[$i])==0 || $value!=0) #ewald cannot have a prime factor other than 2,3,5 and has to be even
 {
  $num++; #oldnumber + 1 = new number
  $i=0; #new number so $i has to become 0 so all the earlier  prime numbers can be checked
 }
 else
 {
  $i++;
 }
}
open(HIGHNUM,">$upload_dir/$filename-highnum.str");
print HIGHNUM <<HIGH_NUM;
*Ewald value
*

set ewaldval $num

return
HIGH_NUM
