/* This functions are to be used with make_report.pl script for SampleFingerprinting workflow @OICR*/


function countDonors(chks,ops,limit){
 var checked = 0
 var donors = []
 for (i=0; i<chks; i++) {
  if (document.getElementById("check" + i).checked) {
    if (checked >= limit ) {
     alert("Too many donors checked, please choose no more than " + limit)
     return
    }
    checked++
    donors.push(document.getElementById("check" + i).value)
  }
 }
 var command = ""
 for (d=0; d<donors.length; d++) {
   d == 0 ? command = command.concat(donors[d]) : command = command.concat("," + donors[d]) 
 }
 document.genrep.donors.value = command
}

