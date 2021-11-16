
<style type="text/css">/********** GMRI Rmarkdown Core Style Sheet - Do Not Modify!!! **********/


/********** Begin Style Sheet **********/

/* Avenir Font from Fonts.com for GMRI Branding */

@import url("http://fast.fonts.net/t/1.css?apiType=css&projectid=806f61f6-d695-4965-a878-820b50bc0269");
    @font-face{
        font-family:"Avenir LT W01_35 Light1475496";
        src:url("Fonts/0078f486-8e52-42c0-ad81-3c8d3d43f48e.woff2") format("woff2"),url("Fonts/908c4810-64db-4b46-bb8e-823eb41f68c0.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_35 Light_1475502";
        src:url("Fonts/a59168c1-917d-4de9-a244-0316c057c357.woff2") format("woff2"),url("Fonts/6dc0e7d8-9284-44e1-8f05-984a41daa3a4.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_45 Book1475508";
        src:url("Fonts/065a6b14-b2cc-446e-9428-271c570df0d9.woff2") format("woff2"),url("Fonts/65d75eb0-2601-4da5-a9a4-9ee67a470a59.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_45 Book O1475514";
        src:url("Fonts/476612d9-282d-4f76-95cd-b4dd31e7ed21.woff2") format("woff2"),url("Fonts/f1ebae2b-5296-4244-8771-5f40e60a564a.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_55 Roman1475520";
        src:url("Fonts/b290e775-e0f9-4980-914b-a4c32a5e3e36.woff2") format("woff2"),url("Fonts/4b978f72-bb48-46c3-909a-2a8cd2f8819c.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_55 Obliqu1475526";
        src:url("Fonts/1a7173fa-062b-49ad-9915-bc57d3bfc1f5.woff2") format("woff2"),url("Fonts/cdda031e-26e9-4269-83d1-5a218caa10db.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_65 Medium1475532";
        src:url("Fonts/17b90ef5-b63f-457b-a981-503bb7afe3c0.woff2") format("woff2"),url("Fonts/c9aeeabd-dd65-491d-b4be-3e0db9ae47a0.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_65 Medium1475538";
        src:url("Fonts/deb5e718-7abb-4df3-9365-edfa95317090.woff2") format("woff2"),url("Fonts/04801919-17ee-4c6b-8b17-eb1965cb3ed6.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_85 Heavy1475544";
        src:url("Fonts/d513e15e-8f35-4129-ad05-481815e52625.woff2") format("woff2"),url("Fonts/61bd362e-7162-46bd-b67e-28f366c4afbe.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_85 Heavy_1475550";
        src:url("Fonts/3c210c80-960f-4684-850b-25390b4d08af.woff2") format("woff2"),url("Fonts/cb5c71ad-e582-4d00-929c-67fbfaeb1c27.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_95 Black1475556";
        src:url("Fonts/c78eb7af-a1c8-4892-974b-52379646fef4.woff2") format("woff2"),url("Fonts/75b36c58-2a02-4057-a537-09af0832ae46.woff") format("woff");
    }
    @font-face{
        font-family:"Avenir LT W01_95 Black_1475562";
        src:url("Fonts/a2477e08-09d9-4d4b-97a9-23a1e22cb44c.woff2") format("woff2"),url("Fonts/19d12bba-92b1-43ad-9bab-cd36a4195c2a.woff") format("woff");
    }





/* PRE-Avenir Fonts: Lato + Raleway font import from google fonts */
@import url('https://fonts.googleapis.com/css?family=Lato');
@import url('https://fonts.googleapis.com/css?family=Raleway&display=swap');

/* add font families as needed: font-family: 'Lato', sans-serif; */


/* Level 1 Headers */
h1 { text-align: left;
     margin: 10px 0 15px 0; 
     font-size: 38px; 
     font-family: Avenir;
} 


/* Headers 2 - 6 */
h2, h3, h4, h5, h6 { 
    color: #333333; 
    margin: 20px 0 5px 0; 
    text-align: left;
    font-family: Avenir;}


/* Sizing/font For Each Header Type */
h2 { font-size: 24px; }
h3 { font-size: 20px; }
h4 { font-size: 18px; }
h5 { font-size: 16px; font-weight: normal; color: #3069aa; text-decoration: underline;}
h6 { font-size: 14px; font-weight: normal; color: #3069aa; }


/* Paragraph Text */
p, ol { margin-top: 10px; 
    font-family: 'Raleway', sans-serif;}


/* Title Author and Date Headers */
h1.title.toc-ignore {margin-top: 10px;}
h4.author, h4.date {
    color: rgb(0,115,109);
    margin-top: 0px; 
    margin-bottom: 5px; 
    font-size: 12px;}



/* Add Spacing Above Headers */
h1, .h1, h2, .h2, h3, .h3 {
    margin-top: 40px;
}


/* Links */
a {
    color: rgb(234,79,18)
}


/***********************************************/


/********  Table of Contents  **********/

/* Highlighted TOC Element */
.list-group-item.active, .list-group-item.active:focus, .list-group-item.active:hover {
    z-index: 2;
    color: #fff;
    background-color: rgb(0,96,138);
    border-color: rgb(0,96,138);
}

/* Default TOC Elements */
.list-group-item, .list-group-item:focus, .list-group-item:hover {
    z-index: 2;
    color: rgb(0,96,138);
    background-color: #fff;
    border-color: rgb(0,96,138);
}


/********  Tab Panels  **********/

/* Navigation Tabs - Highlighted Tabset Pills */
.nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus {
    color: #fff;
    background-color: rgb(0,115,109) ;
    }

/* Navigation Tabs - Default Tabset Pills */
.nav-pills > li > a, .nav-pills > li > a:hover, .nav-pills > li > a:focus {
    color: rgb(0,115,109);
    background-color: #fff;
    }


/* Second Level Tabs - Active */
.nav.nav-tabs > li.active  > a, .nav-tabs > li.active > a:hover, .nav-tabs > li.active > a:focus {
    color: #fff;
    background-color: rgb(83,83,83) ;
    }

/* Second Level Tabs - Inactive */
.nav.nav-tabs > li  > a, .nav-tabs > li > a:hover, .nav-tabs > li > a:focus {
    color: rgb(83,83,83);
    background-color: #fff;
    }



/********** End Core Style Sheet **********/
/********** End Core Style Sheet **********/



</style>

# NEFSC Trawl Survey Analyses

This repository is for analyses relating to the NEFSC bottom trawl
survey and the Maine & New Hampshire trawl survey. Research topics
relate to size-spectrum analyses and species distribution modeling.

## Repository Contents and Reports

Links to other resources and reports detailing EDA, troubleshooting, and
analysis outcomes for NMFS trawl data can be found at the github pages
[index](https://adamkemberling.github.io/nefsc_trawl/).

## Data Preparation

Data for these analyses was obtained from the Northeast Fisheries
Science Center. Details on how data is cleaned prior to any analyses can
be found in our [standard cleanup
doc](https://adamkemberling.github.io/nefsc_trawl/01_Survdat_Standard_Cleanup.html)

<!DOCTYPE html>
<!---- Rmarkdown Custom Footer for Author Details - Modify personal info as needed ------>

<!---- Custom Footer ------->
&nbsp;

<!---- Add GMRI Logo to footer ---->

<!--
<p style="text-align: center;">
    <img src="presentations/gmri_logo.png" alt="GMRI Logo" height="91" width="248">
</p>
-->

<!---- Author Details  ------>
<hr />
<p style="text-align: center; font-size: 60%; padding: 0px;">A work by <a href="https://github.com/adamkemberling/">Adam A. Kemberling</a></p>
<p style="text-align: center; font-size: 60%; padding: 0px;"><span style="color: #808080;"><em>Akemberling@gmri.org</em></span></p>





<!-- Add Font-Awesome Icon Library -->
<!--

<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css">

-->

<!-- Add Icons with Personal Links -->
<!--

<p style="text-align: center;">
    <a href="https://twitter.com/r_graph_gallery?lang=en" class="fa fa-twitter"></a>
    <a href="https://www.linkedin.com/in/yan-holtz-2477534a/" class="fa fa-linkedin"></a>
    <a href="https://github.com/holtzy/" class="fa fa-github"></a>
</p>

-->

&nbsp;

<!--- End Custom Footer --->

</html>
