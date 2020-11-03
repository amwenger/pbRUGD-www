<?php
  $filename = $_POST["filename"];
  $filedata = $_POST["filedata"];
  $filedata = str_replace("<EOL>", PHP_EOL, $filedata);

  header("Content-type: text/plain; charset=UTF-8");
  header("Content-Disposition: attachment; filename=" . $filename);
  header("Expires: 0");
  header("Cache-Control: must-revalidate");
  header("Pragma: public");
  header("Content-Length: " . strval(strlen($filedata)));
  print($filedata);
?>
