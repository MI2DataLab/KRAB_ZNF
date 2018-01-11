-   [Intro](#intro)
-   [Data download and preparation](#data-download-and-preparation)
-   [Data analysis](#data-analysis)
    -   [Categorisation of gene
        expressions](#categorisation-of-gene-expressions)
    -   [Log-rank tests](#log-rank-tests)
    -   [Examples of Kaplan-Meier survival curves for some interesting
        cases](#examples-of-kaplan-meier-survival-curves-for-some-interesting-cases)
        -   [GBMLGG/ZNF485/median](#gbmlggznf485median)
        -   [KIRC/ZNF468/median](#kircznf468median)

Intro
=====

This report provides an overview of survival analysis of cancer patients
based on the TCGA data. The aim was to identify pairs of (cancer cohort,
gene) which possibly indicate the difference in survival time in groups
of patients with high and low expression of a given gene.

TCGA data used for the analysis is obtained via RTCGA family of R
packages.


[Link to Shiny app hosted on
shinyapps.io](https://rafalcyl.shinyapps.io/krab_survival_analysis/) \#

Data download and preparation
=============================


Packages used:

    library(RTCGA)          # Functions for querying the data 
    library(RTCGA.rnaseq)   # Gene expression data
    library(RTCGA.clinical) # Clinical data
    library(survival)       # Survival analysis tools
    library(survminer)      # Survival analysis tools
    library(tidyverse)
    library(knitr)
    library(pheatmap)

Steps of data preparation (R script avialable in /scripts)

1.  Genenes selected for investigation

2.  Extractiong column names for selected genes
3.  Downloading expressions
4.  Downloading survival times
5.  Joining survival times with gene expressions
6.  Check for join multiplicates
7.  We remove SKCM.rnaseq, because there is only one observation.
8.  Shorten gene names
9.  Data overview

<table>
<thead>
<tr class="header">
<th align="left">dataset</th>
<th align="right">n_of_observations</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">ACC</td>
<td align="right">79</td>
</tr>
<tr class="even">
<td align="left">BLCA</td>
<td align="right">403</td>
</tr>
<tr class="odd">
<td align="left">BRCA</td>
<td align="right">1091</td>
</tr>
<tr class="even">
<td align="left">CESC</td>
<td align="right">304</td>
</tr>
<tr class="odd">
<td align="left">CHOL</td>
<td align="right">36</td>
</tr>
<tr class="even">
<td align="left">COAD</td>
<td align="right">282</td>
</tr>
<tr class="odd">
<td align="left">COADREAD</td>
<td align="right">376</td>
</tr>
<tr class="even">
<td align="left">DLBC</td>
<td align="right">15</td>
</tr>
<tr class="odd">
<td align="left">ESCA</td>
<td align="right">183</td>
</tr>
<tr class="even">
<td align="left">GBM</td>
<td align="right">152</td>
</tr>
<tr class="odd">
<td align="left">GBMLGG</td>
<td align="right">663</td>
</tr>
<tr class="even">
<td align="left">HNSC</td>
<td align="right">518</td>
</tr>
<tr class="odd">
<td align="left">KICH</td>
<td align="right">64</td>
</tr>
<tr class="even">
<td align="left">KIPAN</td>
<td align="right">883</td>
</tr>
<tr class="odd">
<td align="left">KIRC</td>
<td align="right">532</td>
</tr>
<tr class="even">
<td align="left">KIRP</td>
<td align="right">287</td>
</tr>
<tr class="odd">
<td align="left">LGG</td>
<td align="right">511</td>
</tr>
<tr class="even">
<td align="left">LIHC</td>
<td align="right">369</td>
</tr>
<tr class="odd">
<td align="left">LUAD</td>
<td align="right">496</td>
</tr>
<tr class="even">
<td align="left">LUSC</td>
<td align="right">487</td>
</tr>
<tr class="odd">
<td align="left">OV</td>
<td align="right">300</td>
</tr>
<tr class="even">
<td align="left">PAAD</td>
<td align="right">178</td>
</tr>
<tr class="odd">
<td align="left">PCPG</td>
<td align="right">179</td>
</tr>
<tr class="even">
<td align="left">PRAD</td>
<td align="right">497</td>
</tr>
<tr class="odd">
<td align="left">READ</td>
<td align="right">94</td>
</tr>
<tr class="even">
<td align="left">SARC</td>
<td align="right">257</td>
</tr>
<tr class="odd">
<td align="left">STAD</td>
<td align="right">36</td>
</tr>
<tr class="even">
<td align="left">STES</td>
<td align="right">219</td>
</tr>
<tr class="odd">
<td align="left">TGCT</td>
<td align="right">133</td>
</tr>
<tr class="even">
<td align="left">THCA</td>
<td align="right">500</td>
</tr>
<tr class="odd">
<td align="left">THYM</td>
<td align="right">119</td>
</tr>
<tr class="even">
<td align="left">UCEC</td>
<td align="right">174</td>
</tr>
<tr class="odd">
<td align="left">UCS</td>
<td align="right">57</td>
</tr>
<tr class="even">
<td align="left">UVM</td>
<td align="right">80</td>
</tr>
</tbody>
</table>

Data analysis
=============

Categorisation of gene expressions
----------------------------------

In order to use log-rank test for comparing two survival curves we need
to split gene expressions variables into categorical variables (high/low
expression).

For our analysis we use two methods for splitting. One is to use median
as a cutpoint and second is to use a cutpoint selected using maximally
selected rank statistics [maximally selected rank
statistics](https://cran.r-project.org/web/packages/maxstat/vignettes/maxstat.pdf)

Note that each cutpoint is calculated in groups (for each combination of
gene/cohort).

Cutpoints obtained using maximally selected rank statistics:

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">ZNF695</th>
<th align="right">ZNF643</th>
<th align="right">ZNF789</th>
<th align="right">ZNF320</th>
<th align="right">ZNF273</th>
<th align="right">ZNF707</th>
<th align="right">ZNF205</th>
<th align="right">ZNF468</th>
<th align="right">ZNF714</th>
<th align="right">ZNF485</th>
<th align="right">ZNF525</th>
<th align="right">ZNF267</th>
<th align="right">ZNF282</th>
<th align="right">ZNF114</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>ACC</td>
<td align="right">13.5524</td>
<td align="right">17.1038</td>
<td align="right">138.9831</td>
<td align="right">479.7433</td>
<td align="right">64.9918</td>
<td align="right">405.7300</td>
<td align="right">528.1063</td>
<td align="right">186.7443</td>
<td align="right">128.3237</td>
<td align="right">34.9286</td>
<td align="right">125.9982</td>
<td align="right">149.7849</td>
<td align="right">777.2845</td>
<td align="right">3.8241</td>
</tr>
<tr class="even">
<td>BLCA</td>
<td align="right">37.6285</td>
<td align="right">25.1572</td>
<td align="right">72.5894</td>
<td align="right">1134.3720</td>
<td align="right">130.4807</td>
<td align="right">205.2003</td>
<td align="right">208.2123</td>
<td align="right">221.7537</td>
<td align="right">105.1213</td>
<td align="right">100.3836</td>
<td align="right">186.4669</td>
<td align="right">439.9510</td>
<td align="right">1294.8975</td>
<td align="right">60.8650</td>
</tr>
<tr class="odd">
<td>BRCA</td>
<td align="right">41.1877</td>
<td align="right">32.0919</td>
<td align="right">110.4221</td>
<td align="right">679.5941</td>
<td align="right">208.7387</td>
<td align="right">161.7473</td>
<td align="right">304.6032</td>
<td align="right">710.9228</td>
<td align="right">396.8301</td>
<td align="right">152.2222</td>
<td align="right">337.4101</td>
<td align="right">248.1545</td>
<td align="right">501.1241</td>
<td align="right">24.7530</td>
</tr>
<tr class="even">
<td>CESC</td>
<td align="right">133.2470</td>
<td align="right">56.0811</td>
<td align="right">137.1003</td>
<td align="right">354.1958</td>
<td align="right">142.6365</td>
<td align="right">197.4919</td>
<td align="right">247.2986</td>
<td align="right">332.9640</td>
<td align="right">108.3196</td>
<td align="right">33.7171</td>
<td align="right">222.4658</td>
<td align="right">485.8830</td>
<td align="right">969.8737</td>
<td align="right">69.7407</td>
</tr>
<tr class="odd">
<td>CHOL</td>
<td align="right">5.7997</td>
<td align="right">9.6878</td>
<td align="right">61.9868</td>
<td align="right">308.3155</td>
<td align="right">80.7099</td>
<td align="right">152.9635</td>
<td align="right">473.0303</td>
<td align="right">351.0577</td>
<td align="right">91.2126</td>
<td align="right">96.3205</td>
<td align="right">152.6598</td>
<td align="right">323.0950</td>
<td align="right">776.6323</td>
<td align="right">1.0488</td>
</tr>
<tr class="even">
<td>COAD</td>
<td align="right">93.7091</td>
<td align="right">77.6899</td>
<td align="right">82.1918</td>
<td align="right">638.4810</td>
<td align="right">130.8467</td>
<td align="right">246.3273</td>
<td align="right">226.6601</td>
<td align="right">527.7651</td>
<td align="right">274.5755</td>
<td align="right">48.7546</td>
<td align="right">214.5125</td>
<td align="right">366.7762</td>
<td align="right">1287.5949</td>
<td align="right">0.8839</td>
</tr>
<tr class="odd">
<td>COADREAD</td>
<td align="right">93.7091</td>
<td align="right">39.7709</td>
<td align="right">82.1918</td>
<td align="right">638.4810</td>
<td align="right">130.8467</td>
<td align="right">246.6960</td>
<td align="right">314.7574</td>
<td align="right">553.3118</td>
<td align="right">274.9208</td>
<td align="right">114.1414</td>
<td align="right">214.5125</td>
<td align="right">360.1799</td>
<td align="right">1257.8526</td>
<td align="right">30.6074</td>
</tr>
<tr class="even">
<td>READ</td>
<td align="right">117.7820</td>
<td align="right">81.9518</td>
<td align="right">174.5574</td>
<td align="right">255.0505</td>
<td align="right">164.4089</td>
<td align="right">320.6215</td>
<td align="right">276.6005</td>
<td align="right">534.3761</td>
<td align="right">94.1638</td>
<td align="right">91.4318</td>
<td align="right">276.1299</td>
<td align="right">266.6256</td>
<td align="right">1034.9498</td>
<td align="right">0.8258</td>
</tr>
<tr class="odd">
<td>DLBC</td>
<td align="right">21.7647</td>
<td align="right">34.2703</td>
<td align="right">82.7878</td>
<td align="right">565.4618</td>
<td align="right">273.2342</td>
<td align="right">251.4440</td>
<td align="right">36.2601</td>
<td align="right">408.4911</td>
<td align="right">214.8481</td>
<td align="right">96.6543</td>
<td align="right">371.3057</td>
<td align="right">1119.5204</td>
<td align="right">676.9349</td>
<td align="right">0.9180</td>
</tr>
<tr class="even">
<td>ESCA</td>
<td align="right">43.2629</td>
<td align="right">87.7386</td>
<td align="right">86.6163</td>
<td align="right">122.2640</td>
<td align="right">167.4192</td>
<td align="right">280.6892</td>
<td align="right">308.9942</td>
<td align="right">629.2933</td>
<td align="right">509.1206</td>
<td align="right">130.4556</td>
<td align="right">108.7304</td>
<td align="right">288.0734</td>
<td align="right">745.9119</td>
<td align="right">6.6735</td>
</tr>
<tr class="odd">
<td>STES</td>
<td align="right">43.2629</td>
<td align="right">87.7386</td>
<td align="right">60.1802</td>
<td align="right">392.6820</td>
<td align="right">167.4192</td>
<td align="right">280.6892</td>
<td align="right">308.9942</td>
<td align="right">725.2935</td>
<td align="right">509.1206</td>
<td align="right">134.0046</td>
<td align="right">314.8561</td>
<td align="right">602.4390</td>
<td align="right">758.6085</td>
<td align="right">0.9117</td>
</tr>
<tr class="even">
<td>LUAD</td>
<td align="right">1.8570</td>
<td align="right">24.9629</td>
<td align="right">37.0882</td>
<td align="right">710.4208</td>
<td align="right">116.6395</td>
<td align="right">168.0726</td>
<td align="right">455.2404</td>
<td align="right">271.1461</td>
<td align="right">157.7564</td>
<td align="right">35.8901</td>
<td align="right">275.9114</td>
<td align="right">248.1969</td>
<td align="right">838.5224</td>
<td align="right">17.7906</td>
</tr>
<tr class="odd">
<td>KIPAN</td>
<td align="right">8.8872</td>
<td align="right">29.7063</td>
<td align="right">101.0249</td>
<td align="right">630.9891</td>
<td align="right">63.6126</td>
<td align="right">145.6241</td>
<td align="right">255.2322</td>
<td align="right">526.6913</td>
<td align="right">44.9882</td>
<td align="right">54.0719</td>
<td align="right">331.7422</td>
<td align="right">278.5700</td>
<td align="right">693.8868</td>
<td align="right">31.4334</td>
</tr>
<tr class="even">
<td>KIRC</td>
<td align="right">7.2533</td>
<td align="right">29.8953</td>
<td align="right">100.7260</td>
<td align="right">650.1506</td>
<td align="right">83.5884</td>
<td align="right">151.3525</td>
<td align="right">414.2466</td>
<td align="right">284.6163</td>
<td align="right">47.7787</td>
<td align="right">56.1694</td>
<td align="right">162.1303</td>
<td align="right">399.1361</td>
<td align="right">813.7438</td>
<td align="right">31.4334</td>
</tr>
<tr class="odd">
<td>GBM</td>
<td align="right">24.8343</td>
<td align="right">37.8021</td>
<td align="right">147.6461</td>
<td align="right">721.1312</td>
<td align="right">215.0033</td>
<td align="right">105.3985</td>
<td align="right">361.1968</td>
<td align="right">274.0538</td>
<td align="right">108.0711</td>
<td align="right">48.6155</td>
<td align="right">77.6119</td>
<td align="right">324.4755</td>
<td align="right">680.3279</td>
<td align="right">25.9230</td>
</tr>
<tr class="even">
<td>GBMLGG</td>
<td align="right">3.4602</td>
<td align="right">32.7771</td>
<td align="right">169.8961</td>
<td align="right">450.5983</td>
<td align="right">277.4086</td>
<td align="right">184.0187</td>
<td align="right">535.9316</td>
<td align="right">192.6821</td>
<td align="right">249.0842</td>
<td align="right">31.7726</td>
<td align="right">108.5098</td>
<td align="right">287.7514</td>
<td align="right">889.4472</td>
<td align="right">99.3007</td>
</tr>
<tr class="odd">
<td>LGG</td>
<td align="right">2.3522</td>
<td align="right">84.7917</td>
<td align="right">74.4585</td>
<td align="right">456.6511</td>
<td align="right">277.4086</td>
<td align="right">231.3261</td>
<td align="right">554.2715</td>
<td align="right">196.1901</td>
<td align="right">364.5437</td>
<td align="right">31.7187</td>
<td align="right">106.3906</td>
<td align="right">300.4222</td>
<td align="right">745.9939</td>
<td align="right">44.3367</td>
</tr>
<tr class="even">
<td>HNSC</td>
<td align="right">33.6451</td>
<td align="right">28.1343</td>
<td align="right">37.5370</td>
<td align="right">56.7792</td>
<td align="right">47.8319</td>
<td align="right">494.0919</td>
<td align="right">524.3965</td>
<td align="right">240.6493</td>
<td align="right">194.8553</td>
<td align="right">78.7037</td>
<td align="right">25.1667</td>
<td align="right">327.9003</td>
<td align="right">1427.3537</td>
<td align="right">35.6537</td>
</tr>
<tr class="odd">
<td>KICH</td>
<td align="right">13.2190</td>
<td align="right">3.5920</td>
<td align="right">58.1717</td>
<td align="right">574.7863</td>
<td align="right">143.0905</td>
<td align="right">159.7836</td>
<td align="right">200.4667</td>
<td align="right">854.1687</td>
<td align="right">407.7634</td>
<td align="right">10.5226</td>
<td align="right">636.8980</td>
<td align="right">114.6461</td>
<td align="right">538.0747</td>
<td align="right">1.5760</td>
</tr>
<tr class="even">
<td>KIRP</td>
<td align="right">12.9862</td>
<td align="right">11.5931</td>
<td align="right">51.4706</td>
<td align="right">434.9385</td>
<td align="right">45.4223</td>
<td align="right">245.3852</td>
<td align="right">426.2993</td>
<td align="right">249.4452</td>
<td align="right">28.6134</td>
<td align="right">53.8533</td>
<td align="right">101.0962</td>
<td align="right">278.5700</td>
<td align="right">1253.5298</td>
<td align="right">0.8111</td>
</tr>
<tr class="odd">
<td>LIHC</td>
<td align="right">0.5301</td>
<td align="right">37.2534</td>
<td align="right">38.0854</td>
<td align="right">281.2555</td>
<td align="right">50.7890</td>
<td align="right">184.0722</td>
<td align="right">261.5912</td>
<td align="right">240.4130</td>
<td align="right">57.1561</td>
<td align="right">20.8535</td>
<td align="right">125.9669</td>
<td align="right">81.1787</td>
<td align="right">701.1446</td>
<td align="right">1.4556</td>
</tr>
<tr class="even">
<td>LUSC</td>
<td align="right">69.3284</td>
<td align="right">23.6132</td>
<td align="right">41.0610</td>
<td align="right">286.3764</td>
<td align="right">174.3970</td>
<td align="right">139.9980</td>
<td align="right">206.8892</td>
<td align="right">309.2764</td>
<td align="right">292.7715</td>
<td align="right">79.0375</td>
<td align="right">387.7286</td>
<td align="right">595.5176</td>
<td align="right">568.6803</td>
<td align="right">2.6326</td>
</tr>
<tr class="odd">
<td>OV</td>
<td align="right">184.3391</td>
<td align="right">19.6780</td>
<td align="right">159.2112</td>
<td align="right">753.3373</td>
<td align="right">184.8124</td>
<td align="right">181.5728</td>
<td align="right">299.6024</td>
<td align="right">330.3611</td>
<td align="right">217.9786</td>
<td align="right">109.1875</td>
<td align="right">242.5349</td>
<td align="right">303.0141</td>
<td align="right">1610.9966</td>
<td align="right">2.2526</td>
</tr>
<tr class="even">
<td>PAAD</td>
<td align="right">0.8180</td>
<td align="right">16.8350</td>
<td align="right">133.5321</td>
<td align="right">470.4218</td>
<td align="right">59.5122</td>
<td align="right">264.1904</td>
<td align="right">161.3679</td>
<td align="right">532.9396</td>
<td align="right">60.4555</td>
<td align="right">30.8552</td>
<td align="right">323.5159</td>
<td align="right">274.9100</td>
<td align="right">569.0578</td>
<td align="right">16.3432</td>
</tr>
<tr class="odd">
<td>PCPG</td>
<td align="right">4.0733</td>
<td align="right">32.7511</td>
<td align="right">74.9263</td>
<td align="right">475.3684</td>
<td align="right">119.4770</td>
<td align="right">141.6630</td>
<td align="right">231.1759</td>
<td align="right">133.5846</td>
<td align="right">208.8963</td>
<td align="right">39.3500</td>
<td align="right">131.9561</td>
<td align="right">199.1814</td>
<td align="right">779.9305</td>
<td align="right">2.9142</td>
</tr>
<tr class="even">
<td>PRAD</td>
<td align="right">22.4776</td>
<td align="right">19.5080</td>
<td align="right">182.6134</td>
<td align="right">576.0300</td>
<td align="right">66.1600</td>
<td align="right">213.9535</td>
<td align="right">593.3990</td>
<td align="right">422.6121</td>
<td align="right">213.6247</td>
<td align="right">75.3473</td>
<td align="right">299.3777</td>
<td align="right">153.1117</td>
<td align="right">929.8311</td>
<td align="right">5.7899</td>
</tr>
<tr class="odd">
<td>SARC</td>
<td align="right">25.5573</td>
<td align="right">39.3736</td>
<td align="right">120.5788</td>
<td align="right">291.9021</td>
<td align="right">96.7514</td>
<td align="right">311.2401</td>
<td align="right">442.5071</td>
<td align="right">261.7870</td>
<td align="right">238.1857</td>
<td align="right">30.7626</td>
<td align="right">124.3842</td>
<td align="right">534.3639</td>
<td align="right">1142.2187</td>
<td align="right">31.4039</td>
</tr>
<tr class="even">
<td>STAD</td>
<td align="right">4.5858</td>
<td align="right">67.9814</td>
<td align="right">91.6351</td>
<td align="right">350.4587</td>
<td align="right">112.4324</td>
<td align="right">208.8332</td>
<td align="right">225.4614</td>
<td align="right">587.2576</td>
<td align="right">421.7437</td>
<td align="right">47.3870</td>
<td align="right">90.5808</td>
<td align="right">285.9659</td>
<td align="right">758.6085</td>
<td align="right">0.7024</td>
</tr>
<tr class="odd">
<td>TGCT</td>
<td align="right">240.9905</td>
<td align="right">224.3893</td>
<td align="right">169.7661</td>
<td align="right">776.5858</td>
<td align="right">687.9664</td>
<td align="right">531.8011</td>
<td align="right">418.7966</td>
<td align="right">612.3937</td>
<td align="right">529.3843</td>
<td align="right">71.8284</td>
<td align="right">955.4475</td>
<td align="right">523.7698</td>
<td align="right">1094.9432</td>
<td align="right">69.8529</td>
</tr>
<tr class="even">
<td>THCA</td>
<td align="right">3.0892</td>
<td align="right">43.9018</td>
<td align="right">48.2166</td>
<td align="right">543.6090</td>
<td align="right">129.5294</td>
<td align="right">160.8240</td>
<td align="right">270.9429</td>
<td align="right">235.9975</td>
<td align="right">391.8152</td>
<td align="right">22.6555</td>
<td align="right">100.4098</td>
<td align="right">396.0519</td>
<td align="right">685.9160</td>
<td align="right">203.0899</td>
</tr>
<tr class="odd">
<td>THYM</td>
<td align="right">81.2516</td>
<td align="right">18.2904</td>
<td align="right">110.1157</td>
<td align="right">278.2516</td>
<td align="right">159.7611</td>
<td align="right">297.3554</td>
<td align="right">439.3728</td>
<td align="right">276.1916</td>
<td align="right">44.0478</td>
<td align="right">54.3950</td>
<td align="right">254.7183</td>
<td align="right">421.5184</td>
<td align="right">865.4344</td>
<td align="right">10.2301</td>
</tr>
<tr class="even">
<td>UCEC</td>
<td align="right">119.1599</td>
<td align="right">21.0000</td>
<td align="right">176.1566</td>
<td align="right">369.6185</td>
<td align="right">75.5562</td>
<td align="right">244.8759</td>
<td align="right">333.9596</td>
<td align="right">292.4394</td>
<td align="right">372.3449</td>
<td align="right">36.9928</td>
<td align="right">398.6388</td>
<td align="right">205.5800</td>
<td align="right">1487.0629</td>
<td align="right">4.0323</td>
</tr>
<tr class="odd">
<td>UCS</td>
<td align="right">253.1902</td>
<td align="right">93.5076</td>
<td align="right">139.8797</td>
<td align="right">378.2074</td>
<td align="right">219.3362</td>
<td align="right">270.7260</td>
<td align="right">446.4434</td>
<td align="right">235.7111</td>
<td align="right">115.0885</td>
<td align="right">40.0330</td>
<td align="right">100.2342</td>
<td align="right">172.2320</td>
<td align="right">655.7462</td>
<td align="right">13.4012</td>
</tr>
<tr class="even">
<td>UVM</td>
<td align="right">9.9267</td>
<td align="right">12.0539</td>
<td align="right">78.7935</td>
<td align="right">411.4088</td>
<td align="right">99.3870</td>
<td align="right">641.4304</td>
<td align="right">324.5295</td>
<td align="right">79.9916</td>
<td align="right">45.6130</td>
<td align="right">68.3012</td>
<td align="right">179.0885</td>
<td align="right">178.8994</td>
<td align="right">1935.8354</td>
<td align="right">47.1850</td>
</tr>
</tbody>
</table>

**Cutpoints obtained using medians**

<table>
<thead>
<tr class="header">
<th align="left">dataset</th>
<th align="right">ZNF695</th>
<th align="right">ZNF643</th>
<th align="right">ZNF789</th>
<th align="right">ZNF320</th>
<th align="right">ZNF273</th>
<th align="right">ZNF707</th>
<th align="right">ZNF205</th>
<th align="right">ZNF468</th>
<th align="right">ZNF714</th>
<th align="right">ZNF485</th>
<th align="right">ZNF525</th>
<th align="right">ZNF267</th>
<th align="right">ZNF282</th>
<th align="right">ZNF114</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">ACC</td>
<td align="right">4.57670</td>
<td align="right">18.26890</td>
<td align="right">65.80890</td>
<td align="right">371.5481</td>
<td align="right">94.87080</td>
<td align="right">260.7744</td>
<td align="right">384.0633</td>
<td align="right">286.60700</td>
<td align="right">144.23450</td>
<td align="right">39.02150</td>
<td align="right">230.6226</td>
<td align="right">161.07380</td>
<td align="right">718.6441</td>
<td align="right">18.38090</td>
</tr>
<tr class="even">
<td align="left">BLCA</td>
<td align="right">55.39690</td>
<td align="right">41.00570</td>
<td align="right">92.35130</td>
<td align="right">611.8058</td>
<td align="right">166.76590</td>
<td align="right">298.0080</td>
<td align="right">340.5401</td>
<td align="right">466.06510</td>
<td align="right">145.90960</td>
<td align="right">60.69800</td>
<td align="right">273.6099</td>
<td align="right">327.49080</td>
<td align="right">1151.1725</td>
<td align="right">17.79700</td>
</tr>
<tr class="odd">
<td align="left">BRCA</td>
<td align="right">35.81420</td>
<td align="right">47.88680</td>
<td align="right">69.76740</td>
<td align="right">537.5079</td>
<td align="right">172.85480</td>
<td align="right">266.2315</td>
<td align="right">304.0768</td>
<td align="right">500.71740</td>
<td align="right">177.85750</td>
<td align="right">86.36190</td>
<td align="right">326.8409</td>
<td align="right">379.66640</td>
<td align="right">728.1288</td>
<td align="right">5.52690</td>
</tr>
<tr class="even">
<td align="left">CESC</td>
<td align="right">66.03905</td>
<td align="right">30.85615</td>
<td align="right">123.86450</td>
<td align="right">212.4275</td>
<td align="right">120.68565</td>
<td align="right">288.1254</td>
<td align="right">396.9140</td>
<td align="right">311.24245</td>
<td align="right">149.00380</td>
<td align="right">44.01705</td>
<td align="right">216.0055</td>
<td align="right">382.92220</td>
<td align="right">927.2707</td>
<td align="right">46.00635</td>
</tr>
<tr class="odd">
<td align="left">CHOL</td>
<td align="right">9.38595</td>
<td align="right">26.83830</td>
<td align="right">95.76355</td>
<td align="right">546.5112</td>
<td align="right">81.34505</td>
<td align="right">279.1318</td>
<td align="right">438.7728</td>
<td align="right">480.88505</td>
<td align="right">86.13525</td>
<td align="right">75.84930</td>
<td align="right">206.3800</td>
<td align="right">299.76035</td>
<td align="right">975.4960</td>
<td align="right">2.19900</td>
</tr>
<tr class="even">
<td align="left">COAD</td>
<td align="right">56.11285</td>
<td align="right">47.81335</td>
<td align="right">112.97765</td>
<td align="right">484.0390</td>
<td align="right">111.53325</td>
<td align="right">227.0582</td>
<td align="right">323.9370</td>
<td align="right">502.32665</td>
<td align="right">153.65200</td>
<td align="right">85.23895</td>
<td align="right">207.4518</td>
<td align="right">350.37955</td>
<td align="right">939.5089</td>
<td align="right">4.07195</td>
</tr>
<tr class="odd">
<td align="left">COADREAD</td>
<td align="right">55.21890</td>
<td align="right">48.91260</td>
<td align="right">114.03790</td>
<td align="right">484.8668</td>
<td align="right">113.82945</td>
<td align="right">231.4810</td>
<td align="right">319.6178</td>
<td align="right">507.81270</td>
<td align="right">152.09125</td>
<td align="right">85.75415</td>
<td align="right">227.9812</td>
<td align="right">355.61775</td>
<td align="right">944.2114</td>
<td align="right">4.12340</td>
</tr>
<tr class="even">
<td align="left">DLBC</td>
<td align="right">17.79800</td>
<td align="right">37.09250</td>
<td align="right">175.18590</td>
<td align="right">353.6854</td>
<td align="right">208.53310</td>
<td align="right">347.1429</td>
<td align="right">38.9312</td>
<td align="right">277.14290</td>
<td align="right">179.04860</td>
<td align="right">87.85710</td>
<td align="right">208.1784</td>
<td align="right">957.33460</td>
<td align="right">788.9845</td>
<td align="right">1.10010</td>
</tr>
<tr class="odd">
<td align="left">ESCA</td>
<td align="right">57.95490</td>
<td align="right">48.15760</td>
<td align="right">94.55080</td>
<td align="right">319.1851</td>
<td align="right">167.41920</td>
<td align="right">191.5474</td>
<td align="right">257.1679</td>
<td align="right">340.86240</td>
<td align="right">306.23330</td>
<td align="right">62.66100</td>
<td align="right">156.3591</td>
<td align="right">478.00390</td>
<td align="right">874.2173</td>
<td align="right">22.63720</td>
</tr>
<tr class="even">
<td align="left">GBM</td>
<td align="right">4.13580</td>
<td align="right">39.45155</td>
<td align="right">133.34765</td>
<td align="right">477.7987</td>
<td align="right">237.76675</td>
<td align="right">165.4587</td>
<td align="right">318.5251</td>
<td align="right">236.68175</td>
<td align="right">178.22870</td>
<td align="right">28.57305</td>
<td align="right">102.4320</td>
<td align="right">273.40645</td>
<td align="right">867.9384</td>
<td align="right">64.24785</td>
</tr>
<tr class="odd">
<td align="left">GBMLGG</td>
<td align="right">8.07440</td>
<td align="right">46.94760</td>
<td align="right">101.02050</td>
<td align="right">401.1903</td>
<td align="right">203.27720</td>
<td align="right">213.3132</td>
<td align="right">345.3393</td>
<td align="right">133.69780</td>
<td align="right">252.18080</td>
<td align="right">37.53350</td>
<td align="right">80.3797</td>
<td align="right">226.55760</td>
<td align="right">857.9147</td>
<td align="right">54.11370</td>
</tr>
<tr class="even">
<td align="left">HNSC</td>
<td align="right">19.87440</td>
<td align="right">37.24730</td>
<td align="right">47.24925</td>
<td align="right">140.1728</td>
<td align="right">101.45950</td>
<td align="right">318.3714</td>
<td align="right">307.0445</td>
<td align="right">261.50635</td>
<td align="right">98.14185</td>
<td align="right">43.26070</td>
<td align="right">135.9057</td>
<td align="right">411.91835</td>
<td align="right">1061.3465</td>
<td align="right">65.46595</td>
</tr>
<tr class="odd">
<td align="left">KICH</td>
<td align="right">4.17970</td>
<td align="right">4.71710</td>
<td align="right">36.45485</td>
<td align="right">709.7993</td>
<td align="right">110.81315</td>
<td align="right">201.7873</td>
<td align="right">364.7158</td>
<td align="right">874.79455</td>
<td align="right">295.58070</td>
<td align="right">15.36045</td>
<td align="right">399.1436</td>
<td align="right">176.22935</td>
<td align="right">702.8828</td>
<td align="right">8.15070</td>
</tr>
<tr class="even">
<td align="left">KIPAN</td>
<td align="right">3.35980</td>
<td align="right">21.95000</td>
<td align="right">87.58410</td>
<td align="right">1020.8184</td>
<td align="right">70.64760</td>
<td align="right">152.3361</td>
<td align="right">361.3642</td>
<td align="right">389.26170</td>
<td align="right">75.71450</td>
<td align="right">43.22770</td>
<td align="right">208.3212</td>
<td align="right">362.59260</td>
<td align="right">669.1561</td>
<td align="right">6.06090</td>
</tr>
<tr class="odd">
<td align="left">KIRC</td>
<td align="right">3.43220</td>
<td align="right">25.85310</td>
<td align="right">87.59545</td>
<td align="right">1110.2584</td>
<td align="right">65.20860</td>
<td align="right">129.0370</td>
<td align="right">328.6197</td>
<td align="right">398.12810</td>
<td align="right">78.49730</td>
<td align="right">52.08850</td>
<td align="right">223.2423</td>
<td align="right">426.51210</td>
<td align="right">600.9410</td>
<td align="right">8.91245</td>
</tr>
<tr class="even">
<td align="left">KIRP</td>
<td align="right">3.17960</td>
<td align="right">17.34510</td>
<td align="right">119.17300</td>
<td align="right">980.7461</td>
<td align="right">82.60390</td>
<td align="right">189.4764</td>
<td align="right">451.2606</td>
<td align="right">336.96240</td>
<td align="right">56.39170</td>
<td align="right">32.95850</td>
<td align="right">161.1271</td>
<td align="right">279.13570</td>
<td align="right">848.1288</td>
<td align="right">2.86940</td>
</tr>
<tr class="odd">
<td align="left">LGG</td>
<td align="right">8.69610</td>
<td align="right">49.39030</td>
<td align="right">93.01470</td>
<td align="right">378.0048</td>
<td align="right">197.27670</td>
<td align="right">228.2514</td>
<td align="right">362.7672</td>
<td align="right">116.12820</td>
<td align="right">275.80230</td>
<td align="right">41.26520</td>
<td align="right">76.4636</td>
<td align="right">214.60780</td>
<td align="right">857.3482</td>
<td align="right">51.81860</td>
</tr>
<tr class="even">
<td align="left">LIHC</td>
<td align="right">1.26280</td>
<td align="right">19.81770</td>
<td align="right">72.46380</td>
<td align="right">139.4813</td>
<td align="right">37.02050</td>
<td align="right">200.1358</td>
<td align="right">318.7615</td>
<td align="right">183.39740</td>
<td align="right">21.47240</td>
<td align="right">39.61670</td>
<td align="right">69.0131</td>
<td align="right">139.14080</td>
<td align="right">759.5517</td>
<td align="right">1.84310</td>
</tr>
<tr class="odd">
<td align="left">LUAD</td>
<td align="right">16.80745</td>
<td align="right">38.58460</td>
<td align="right">82.55375</td>
<td align="right">477.5593</td>
<td align="right">131.83670</td>
<td align="right">232.2317</td>
<td align="right">277.9214</td>
<td align="right">434.31055</td>
<td align="right">180.52735</td>
<td align="right">57.95845</td>
<td align="right">190.1596</td>
<td align="right">391.27730</td>
<td align="right">733.8367</td>
<td align="right">15.25460</td>
</tr>
<tr class="even">
<td align="left">LUSC</td>
<td align="right">52.32940</td>
<td align="right">48.95610</td>
<td align="right">86.57430</td>
<td align="right">287.5008</td>
<td align="right">183.20940</td>
<td align="right">261.3042</td>
<td align="right">246.1825</td>
<td align="right">365.16650</td>
<td align="right">166.26010</td>
<td align="right">51.87820</td>
<td align="right">238.1848</td>
<td align="right">426.17830</td>
<td align="right">931.7450</td>
<td align="right">24.47910</td>
</tr>
<tr class="odd">
<td align="left">OV</td>
<td align="right">82.50060</td>
<td align="right">45.98915</td>
<td align="right">119.56695</td>
<td align="right">457.2094</td>
<td align="right">209.28670</td>
<td align="right">347.8374</td>
<td align="right">446.6301</td>
<td align="right">267.43470</td>
<td align="right">329.97195</td>
<td align="right">91.77680</td>
<td align="right">170.5825</td>
<td align="right">258.14105</td>
<td align="right">1339.3078</td>
<td align="right">10.84365</td>
</tr>
<tr class="even">
<td align="left">PAAD</td>
<td align="right">6.03660</td>
<td align="right">25.96830</td>
<td align="right">103.48495</td>
<td align="right">419.2163</td>
<td align="right">91.60095</td>
<td align="right">238.5402</td>
<td align="right">361.4432</td>
<td align="right">347.60110</td>
<td align="right">81.04045</td>
<td align="right">50.79560</td>
<td align="right">199.1255</td>
<td align="right">445.62825</td>
<td align="right">818.4164</td>
<td align="right">6.51780</td>
</tr>
<tr class="odd">
<td align="left">PCPG</td>
<td align="right">1.06300</td>
<td align="right">21.67180</td>
<td align="right">70.37360</td>
<td align="right">338.0391</td>
<td align="right">56.40310</td>
<td align="right">153.1619</td>
<td align="right">378.2247</td>
<td align="right">114.95480</td>
<td align="right">98.42880</td>
<td align="right">25.10040</td>
<td align="right">129.4569</td>
<td align="right">147.63950</td>
<td align="right">540.1372</td>
<td align="right">2.87150</td>
</tr>
<tr class="even">
<td align="left">PRAD</td>
<td align="right">7.87990</td>
<td align="right">30.24030</td>
<td align="right">122.69550</td>
<td align="right">555.9705</td>
<td align="right">87.33620</td>
<td align="right">220.8901</td>
<td align="right">367.7533</td>
<td align="right">420.77750</td>
<td align="right">138.47330</td>
<td align="right">80.49790</td>
<td align="right">311.6432</td>
<td align="right">155.39070</td>
<td align="right">1075.5814</td>
<td align="right">2.34840</td>
</tr>
<tr class="odd">
<td align="left">READ</td>
<td align="right">51.08200</td>
<td align="right">53.81235</td>
<td align="right">115.67790</td>
<td align="right">494.2496</td>
<td align="right">123.07210</td>
<td align="right">249.0764</td>
<td align="right">301.6316</td>
<td align="right">521.46760</td>
<td align="right">144.61385</td>
<td align="right">87.96295</td>
<td align="right">267.4578</td>
<td align="right">381.66415</td>
<td align="right">959.8034</td>
<td align="right">4.29305</td>
</tr>
<tr class="even">
<td align="left">SARC</td>
<td align="right">13.17470</td>
<td align="right">36.82630</td>
<td align="right">73.18410</td>
<td align="right">363.9500</td>
<td align="right">116.47140</td>
<td align="right">234.6613</td>
<td align="right">403.2935</td>
<td align="right">186.07520</td>
<td align="right">101.71680</td>
<td align="right">33.80280</td>
<td align="right">124.1756</td>
<td align="right">298.10180</td>
<td align="right">934.1085</td>
<td align="right">15.69260</td>
</tr>
<tr class="odd">
<td align="left">STAD</td>
<td align="right">50.39940</td>
<td align="right">51.72035</td>
<td align="right">92.01475</td>
<td align="right">398.5734</td>
<td align="right">112.38095</td>
<td align="right">185.1175</td>
<td align="right">257.1717</td>
<td align="right">465.52385</td>
<td align="right">258.44365</td>
<td align="right">61.06590</td>
<td align="right">209.2262</td>
<td align="right">495.87195</td>
<td align="right">852.9116</td>
<td align="right">5.69955</td>
</tr>
<tr class="even">
<td align="left">STES</td>
<td align="right">57.24000</td>
<td align="right">48.79540</td>
<td align="right">94.35640</td>
<td align="right">335.8714</td>
<td align="right">158.99650</td>
<td align="right">190.7514</td>
<td align="right">257.1679</td>
<td align="right">353.98300</td>
<td align="right">298.25520</td>
<td align="right">61.77610</td>
<td align="right">159.6811</td>
<td align="right">482.31990</td>
<td align="right">863.3507</td>
<td align="right">15.32910</td>
</tr>
<tr class="odd">
<td align="left">TGCT</td>
<td align="right">205.46400</td>
<td align="right">102.68490</td>
<td align="right">148.16500</td>
<td align="right">772.3120</td>
<td align="right">841.64860</td>
<td align="right">359.2767</td>
<td align="right">406.6032</td>
<td align="right">467.16790</td>
<td align="right">438.32360</td>
<td align="right">123.50880</td>
<td align="right">554.3294</td>
<td align="right">489.09300</td>
<td align="right">815.4667</td>
<td align="right">132.54110</td>
</tr>
<tr class="even">
<td align="left">THCA</td>
<td align="right">3.48045</td>
<td align="right">34.93840</td>
<td align="right">73.35835</td>
<td align="right">418.6057</td>
<td align="right">170.60855</td>
<td align="right">218.8245</td>
<td align="right">348.3640</td>
<td align="right">385.03125</td>
<td align="right">387.45980</td>
<td align="right">32.27385</td>
<td align="right">147.7005</td>
<td align="right">276.83075</td>
<td align="right">635.6639</td>
<td align="right">174.54300</td>
</tr>
<tr class="odd">
<td align="left">THYM</td>
<td align="right">34.68460</td>
<td align="right">25.23430</td>
<td align="right">169.04000</td>
<td align="right">395.5590</td>
<td align="right">166.04240</td>
<td align="right">351.8105</td>
<td align="right">341.9270</td>
<td align="right">220.99870</td>
<td align="right">143.08130</td>
<td align="right">50.61690</td>
<td align="right">238.9023</td>
<td align="right">279.02620</td>
<td align="right">981.6066</td>
<td align="right">9.38230</td>
</tr>
<tr class="even">
<td align="left">UCEC</td>
<td align="right">122.25655</td>
<td align="right">40.54370</td>
<td align="right">146.23165</td>
<td align="right">700.2405</td>
<td align="right">164.99865</td>
<td align="right">290.5149</td>
<td align="right">544.7861</td>
<td align="right">360.16900</td>
<td align="right">165.30625</td>
<td align="right">69.44235</td>
<td align="right">242.2034</td>
<td align="right">259.57845</td>
<td align="right">1139.3367</td>
<td align="right">25.78640</td>
</tr>
<tr class="odd">
<td align="left">UCS</td>
<td align="right">147.98400</td>
<td align="right">52.26130</td>
<td align="right">128.24880</td>
<td align="right">552.5890</td>
<td align="right">195.96660</td>
<td align="right">396.1112</td>
<td align="right">546.0260</td>
<td align="right">225.76180</td>
<td align="right">188.28050</td>
<td align="right">75.69910</td>
<td align="right">189.0651</td>
<td align="right">239.21430</td>
<td align="right">1049.6300</td>
<td align="right">29.08590</td>
</tr>
<tr class="even">
<td align="left">UVM</td>
<td align="right">5.73825</td>
<td align="right">20.07800</td>
<td align="right">65.03550</td>
<td align="right">242.7931</td>
<td align="right">46.27870</td>
<td align="right">444.8644</td>
<td align="right">583.3197</td>
<td align="right">82.27815</td>
<td align="right">28.86065</td>
<td align="right">53.82010</td>
<td align="right">73.0601</td>
<td align="right">74.25925</td>
<td align="right">1223.1123</td>
<td align="right">18.95840</td>
</tr>
</tbody>
</table>

Log-rank tests
--------------

Now, for each pair \[gene, cohort\] we calculate p-values of log-rank
tests.

Low p-value indicates that there is a significant difference in survival
time for two groups i.e. high/low expression.

Heatmap of p-values for log-rank test using cutpoint via maximally
selected rank statistics.
![](report/survival_analysis_v1_files/figure-markdown_strict/unnamed-chunk-7-1.png)

Heatmap of p-values for log-rank test using median cutpoint.
![](report/survival_analysis_v1_files/figure-markdown_strict/unnamed-chunk-8-1.png)

Examples of Kaplan-Meier survival curves for some interesting cases
-------------------------------------------------------------------

### GBMLGG/ZNF485/median

In these case p &lt; 0.05 and we observe that survival curves differs
i.e. lower expression = longer survival.

![](report/survival_analysis_v1_files/figure-markdown_strict/unnamed-chunk-9-1.png)

### KIRC/ZNF468/median

In these case p &gt; 0.05 and we observe that survival curves are
realively close to each other.

![](report/survival_analysis_v1_files/figure-markdown_strict/unnamed-chunk-10-1.png)
