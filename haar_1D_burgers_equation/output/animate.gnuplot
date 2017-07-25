if (i<10) {
    plot '000'.i.'.dat' with lines
} else {

}
if (i<100&&i>=10) {
    plot '00'.i.'.dat' with lines
} else {

}    
if (i<1000&&i>=100) {
    plot '0'.i.'.dat' with lines
} else {

}
i=i+1
if (i<n) reread
