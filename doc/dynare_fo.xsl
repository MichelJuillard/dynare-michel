<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:fo="http://www.w3.org/1999/XSL/Format"
                xmlns:doc="http://nwalsh.com/xsl/documentation/1.0"
                exclude-result-prefixes="doc"
                version='1.0'>

<!-- It is important to use indent="no" here, otherwise verbatim -->
<!-- environments get broken by indented tags...at least when the -->
<!-- callout extension is used...at least with some processors -->
<xsl:import href="http://docbook.sourceforge.net/release/xsl/current/fo/docbook.xsl"/>
<xsl:import href="dynare.xsl"/>

<xsl:import href="dynare_lists.xsl"/>

<xsl:param name="passivetex.extensions">1</xsl:param>
<xsl:attribute-set name="list.item.spacing">
  <xsl:attribute name="line-height">0em</xsl:attribute>
</xsl:attribute-set>
<xsl:attribute-set name="compact.list.item.spacing">
  <xsl:attribute name="space-before.optimum">0em</xsl:attribute>
  <xsl:attribute name="space-before.minimum">0em</xsl:attribute>
  <xsl:attribute name="space-before.maximum">0em</xsl:attribute>
  <xsl:attribute name="line-height">0em</xsl:attribute>
</xsl:attribute-set>
 

</xsl:stylesheet>
