function Output = AvgPhaseDiff(PhaseMap, Point1, Point2)
Point1 = squeeze(Point1)';
Point2 = squeeze(Point2)';
Output = -mean(squeeze(PhaseMap(Point1(1),Point1(2),:)-PhaseMap(Point2(1),Point2(1),:)));