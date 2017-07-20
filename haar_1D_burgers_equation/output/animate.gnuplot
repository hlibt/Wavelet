if (i<10) {
    plot '000'.i.'.dat' with lines
} else {
    plot '00'.i.'.dat' with lines
} 
i=i+1
if (i<n) reread
