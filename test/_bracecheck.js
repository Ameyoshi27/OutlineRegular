const fs=require('fs');
const text=fs.readFileSync('src/main.cpp','utf8');
let depth=0,line=1;
let inStr=false,inChar=false,inLine=false,inBlock=false;
let lastZeroLine=0;
for(let i=0;i<text.length;i++){
  const c=text[i],n=text[i+1];
  if(c==='\n'){line++;inLine=false;continue;}
  if(inLine)continue;
  if(inBlock){if(c==='*'&&n==='/'){inBlock=false;i++;}continue;}
  if(inStr){if(c==='\'){i++;continue;}if(c==='"'){inStr=false;}continue;}
  if(inChar){if(c==='\'){i++;continue;}if(c==="'"){inChar=false;}continue;}
  if(c==='/'&&n==='/'){inLine=true;continue;}
  if(c==='/'&&n==='*'){inBlock=true;i++;continue;}
  if(c==='"'){inStr=true;continue;}
  if(c==="'"){inChar=true;continue;}
  if(c==='{'){depth++;}
  if(c==='}'){depth--;if(depth===0)lastZeroLine=line;if(depth<0){console.log('NEGATIVE at line',line);process.exit(0);}}
}
console.log('final depth:',depth,'last zero-depth line:',lastZeroLine,'total:',text.split('\n').length);
